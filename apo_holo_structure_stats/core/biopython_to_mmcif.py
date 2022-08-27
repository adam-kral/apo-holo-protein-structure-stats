import logging
from collections import defaultdict
from typing import Dict, DefaultDict, Tuple, NamedTuple, Iterable, List, Set

from Bio.PDB.Residue import Residue

logger = logging.getLogger(__name__)

class ResidueId(NamedTuple):
    """ Allows unique identification of a residue within a PDB structure

     Esp. within Bio.PDB.Structure (otherwise it'd be better with label_seq_id and label_asym_id) """
    label_seq_id: int

    chain_id: str  # is already present in e.g. ChainResidueData, but usually I expect just SetOfResidueData (may be from multiple chains,
    # simplest, ok)

    @classmethod
    def from_bio_residue(cls, residue: Residue, residue_id_mapping: 'BiopythonToMmcifResidueIds.Mapping'):
        return ResidueId(residue_id_mapping.to_label_seq_id(BioResidueId(*residue.id)), residue.get_parent().id)


class BioResidueId(NamedTuple):
    """ Allows unique identification of a residue within a PDB structure

     Esp. within Bio.PDB.Structure (otherwise it'd be better with label_seq_id and label_asym_id) """
    hetero_flag: str
    auth_seq_id: int
    insertion_code: str
    #
    # chain_id: str  # is already present in e.g. ChainResidueData, but usually I expect just SetOfResidueData (may be from multiple chains,
    # # simplest, ok)

    @classmethod
    def from_bio(cls, residue: Residue):
        return cls(*residue.get_id())


class BiopythonToMmcifResidueIds:
    class Mapping:
        def __init__(self, entity_poly_id: str, label_asym_id: str):
            self.entity_poly_id = entity_poly_id
            self.label_asym_id = label_asym_id

            self.label_seq_id__to__bio_pdb: Dict[int, BioResidueId] = {}
            self.bio_pdb__to__label_seq_id: Dict[BioResidueId, int] = {}

        def add_residue(self, label_seq_id: int, biopython_residue_id: BioResidueId):
            self.label_seq_id__to__bio_pdb[label_seq_id] = biopython_residue_id
            self.bio_pdb__to__label_seq_id[biopython_residue_id] = label_seq_id

        def to_bio(self, label_seq_id: int) -> BioResidueId:
            return self.label_seq_id__to__bio_pdb[label_seq_id]

        def to_label_seq_id(self, bio_residue_id: BioResidueId) -> int:
            return self.bio_pdb__to__label_seq_id[bio_residue_id]

        def all_to_bio(self, label_seq_ids: Iterable[int]) -> Iterable[BioResidueId]:
            return map(lambda label_seq_id: self.to_bio(label_seq_id), label_seq_ids)

        def find_label_seq(self, auth_seq_id: int) -> int:
            """ Similar to `to_label_seq_id`, but ignores insertion codes and hetero flags.

            Use only if you have just auth_seq_id (results in a paper). """
            auth_to_label = {bio_residue_id.auth_seq_id: label_seq
                   for bio_residue_id, label_seq in self.bio_pdb__to__label_seq_id.items()}
            return auth_to_label[auth_seq_id]

    Chains = Dict[str, Mapping]
    Models = DefaultDict[int, Chains]

    EntityPolySequences = DefaultDict[str, Dict[int, str]]  # entity_poly_id -> {label_seq_id: three_letter_AA_code,...}

    @staticmethod
    def get_all_entity_seq(mmcif_dict) -> Tuple[EntityPolySequences, Set[int]]:
        """ Returns dict entity_id -> list of monomer 3-letter codes.

        Will return only polymer entities. At the moment entities with sequence microheterogeneity, as in
        https://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Categories/entity_poly_seq.html
        are skipped.
        """
        entity_seqs: BiopythonToMmcifResidueIds.EntityPolySequences = defaultdict(dict)
        entity_ids_with_seq_heterogeneity = set()  # will be deleted (that's current behavior, or could be returned with the seqs. somehow)

        for entity_id, seq_num, monomer_id in zip(mmcif_dict["_entity_poly_seq.entity_id"],
                                                  mmcif_dict["_entity_poly_seq.num"],
                                                  mmcif_dict["_entity_poly_seq.mon_id"]):

            seq_num = int(seq_num)  # only this field of the three is integer

            if seq_num in entity_seqs[entity_id]:
                logger.warning(f'{mmcif_dict["_entry.id"]} Microheterogeneity in sequence. As in '
                                f'https://mmcif.wwpdb.org/dictionaries/mmcif_std.dic/Categories/entity_poly_seq.html '
                                '(Multiple sequences for the structure.) Skipping...')
                entity_ids_with_seq_heterogeneity.add(entity_id)

            entity_seqs[entity_id][seq_num] = monomer_id

        # skipping entities (chains) with a position where multiple amino acids are possible
        # todo ensure that this is not to have one entity for two different chains with entirely different coordinates
        #  (At the moment I only have one example: 5ZA2, and it does not manifest), in any case, this heterogeneity is
        #  rare
        return entity_seqs, entity_ids_with_seq_heterogeneity
    #
    # def get_chain_mappings(self):
    #     for chain_id in self.label_seq_id__to__bio_pdb:
    #         entity_id = self.chain_id__to__entity_poly_id[chain_id]
    #         label_seq__to__biopython = self.label_seq_id__to__bio_pdb[chain_id]
    #         biopython__to__label_seq = self.bio_pdb__to__label_seq_id[chain_id]
    #
    #         yield entity_id, chain_id, self.entity_poly_sequences[entity_id], label_seq__to__biopython, biopython__to__label_seq

    class BiopythonCompatibilityModelIdCounter:
        """ Translates pdbx_PDB_model_num to Biopython's model id (0-based, no gaps) """
        def __init__(self, mmcif_dict):
            try:
                self.model_id_list = [int(n) for n in mmcif_dict["_atom_site.pdbx_PDB_model_num"]]
            except KeyError:
                # No model number column
                self.model_id_list = None

            self.last_model_id = -1  # pdbx_PDB_model_num, usually 1-based, probably with no guarantees
            self.biopython_model_id = -1  # 0-based, sequential (no gaps)

        def get_for_atom(self, atom_serial):
            # biopython's code:
            if self.model_id_list is not None:
                # model column exists; use it
                model_id = self.model_id_list[atom_serial]
                if self.last_model_id != model_id:
                    # if pdbx_PDB_model_num changes, update it and start new model id
                    self.last_model_id = model_id
                    self.biopython_model_id += 1
            else:
                # no explicit model column

                # probably a bug in Biopython. The model would have id -1.
                # But will this branch be ever used?
                #   _atom_site.pdbx_PDB_model_num: Used in current PDB entries: Yes, in about 100.0 % of entries
                #   probably not...
                logger.warning(f'No pdbx_PDB_model_num column in mmcif. Biopython\'s model would be at id -1. This project code only uses model at index 0. Will crash afterwards.')
                pass

            return self.biopython_model_id

    @classmethod
    def create(cls, mmcif_dict) -> Tuple[Models, EntityPolySequences, Set[int]]:
        # structure as Biopython PDB module: Structure (this class) -> Models -> Chains
        models = cls.Models(dict)  # ugly defaultdict initialization (typing doesn't help)

        entity_poly_sequences, entity_ids_with_seq_heterogeneity = cls.get_all_entity_seq(mmcif_dict)

        # mmcif atom list spec: https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v40.dic/Categories/atom_site.html

        chain_id_list = mmcif_dict["_atom_site.auth_asym_id"]  # not required in spec, but is in 100.0 % entries
        label_asym_id_list = mmcif_dict["_atom_site.label_asym_id"]  # not required in spec, but is in 100.0 % entries
        entity_id_list = mmcif_dict["_atom_site.label_entity_id"]
        resname_list = mmcif_dict["_atom_site.label_comp_id"]
        icode_list = mmcif_dict["_atom_site.pdbx_PDB_ins_code"]  # not required in spec, in about 3.2 % of entries (probably means the field is present but with
        group_PDB_list = mmcif_dict["_atom_site.group_PDB"]  # ATOM or HETATM, not required in spec, but is in 100.0 % entries

        # all list are the same length - number of atoms in _atom_site loop
        total_atoms = len(chain_id_list)

        biopython_model_id_counter = cls.BiopythonCompatibilityModelIdCounter(mmcif_dict)

        # this is the key data that we need and is not present in BioPython.PDB.Residue
        label_seq_id_list = mmcif_dict["_atom_site.label_seq_id"]

        try:
            auth_seq_id_list = mmcif_dict["_atom_site.auth_seq_id"]
        except KeyError:
            # this try-except block matches BioPython's MMCIFParser behavior
            auth_seq_id_list = label_seq_id_list

        # two special chars as placeholders in the mmCIF format
        # for item values that cannot be explicitly assigned
        # see: pdbx/mmcif syntax web page
        _unassigned = {".", "?"}

        last_label_seq_id = None

        for i in range(total_atoms):
            # skip entities we didn't parse in get_all_entity_seq (generally non-polymer entities)
            if entity_id_list[i] not in entity_poly_sequences:
                continue

            # # skip if _atom_site.pdbx_PDB_model_num is not the model we want (if not set, everything is in a single BioPython's model)
            # if model_id_list and model_id_list[i] != pdbx_PDB_model_num:
            #     continue

            # skip HETATMs (not part of polypeptide, label_seq_id is undefined, unassigned)
            if group_PDB_list[i] == "HETATM":
                continue

            label_seq_id = label_seq_id_list[i]

            # skip non-polymer entities
            if label_seq_id in _unassigned:
                continue

            # see if this atom begins a new residue
            if label_seq_id == last_label_seq_id:
                continue

            chain_id = chain_id_list[i]

            # prepare Residue attributes as in BioPython
            resname = resname_list[i]
            int_resseq = int(auth_seq_id_list[i])
            icode = icode_list[i]
            if icode in _unassigned:
                icode = " "

            group_PDB = group_PDB_list[i]
            if group_PDB == "HETATM":
                if resname == "HOH" or resname == "WAT":
                    hetatm_flag = "W"
                else:
                    hetatm_flag = "H"  # This might be used ("W" branch probably won't, as we skip non-polymers),
                    # ex. phosphoserine in 5za2 (which would be skipped though due to microheterogeneity)
            else:
                hetatm_flag = " "

            biopython_residue_id = BioResidueId(hetatm_flag, int_resseq, icode)
            model_id = biopython_model_id_counter.get_for_atom(i)

            # construct the mapping
            model = models[model_id]
            try:
                chain_mapping = model[chain_id]
            except KeyError:
                # todo neni zaruceny, ze auth_asym_id ma jedinou entitu, label_asym_id ano, mozna bych to mohl
                #  kontrolovat/zjistit
                # here I assume that an auth_asym_id chain covers the label_sym_id chain (i.e. label_asym_id is same
                # for all the residues)
                # if it weren't true, I should make one mapping for the whole structure, a mapping of tuples
                # (label/auth chain_id, label/auth residue_id) onto each other.
                model[chain_id] = chain_mapping = BiopythonToMmcifResidueIds.Mapping(entity_id_list[i],
                                                                                     label_asym_id_list[i])

            chain_mapping.add_residue(int(label_seq_id), biopython_residue_id)

            last_label_seq_id = label_seq_id

        return models, entity_poly_sequences, entity_ids_with_seq_heterogeneity


            # [done] add building polypeptide sequence ("microheterogeneity" - have to use resname to resolve ambiguities)
            # teď nevim, co dělat s label_alt_id... a taky sekvence neni povinna
            # navic     _entity_poly.pdbx_seq_one_letter_code_can nema moznost heterogenity jako entity_poly_seq, tak si tam jeden zbytek asi vyberou?
            # viděl jsem jenom příklad 5ZA2 (z prikladu, kdy jde o jiny residue), tak nevim, jak jinde. Tam byl zrovna alt. residue s jinym cmpd_id oznacen jako HETATM, takže
            # by to v biopython parseru nevadilo.

            # nejde alignovat ne-python-stringy, takže to umi jenom one-letter-code NE pry jde, když to je list, ale je to o dost pomalejsi

            # CifSeqresIterator maj napsany spatne. 1) chain id je to label_asym_id, tj. to mmcifovsky (nekompatibilni s modulem BioPython.PDB),
            # 2) heterogenita to rozbije, sekvence je pak delší

            # varianta 2 - udelat to jak rikal, s tim mapovanim na uniprot. Musel bych kontrolovat, ze tam nejsou zadny
            # gapy ani mismatche, to by celkem slo. A pak vzdy vzal minimum z těch rozměrů 2 sekv do toho porovnání.
            # Ale když jsou unobserved, nepoznam, jestli jde o gap, nebo ne, pak bych musel věřit tomu namapovaní, i když
            # by mohlo být špatně - mluvěj tam o seq. identity>95? Takze tam muzou byt i gapy, ktery nemuzu detekovat?
            # a mismatche v tech unobserved bych taky nezdetekoval
            #   - ne to mozna jo, kdybych se podival na unp->pdb, ten segment nebo cely, nevim, nikde nepisou, jestli dovolujou gapy
            #  - je to trochu nepruhledny (jen pisou this complex procedure also enables annotation of differences, such as variants, isoforms, modified residues, microheterogeneity or engineered mutations, between the sample sequence and the UniProtKB sequence. )
            #  -  a srovnavani s druhou sekvenci by bylo slozity a moc nedavalo smysl, kdyz muzu jednoduse srovnat sekvenci z experimentu..
            #  - jenze to ma taky problem



