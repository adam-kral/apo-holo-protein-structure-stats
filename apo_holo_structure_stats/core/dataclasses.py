from bisect import bisect_left, bisect
from dataclasses import dataclass
from enum import Enum, auto
from typing import List, Dict, TypeVar, Generic, Iterator, Iterable, Tuple, Any

import numpy as np
from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue

from apo_holo_structure_stats.core.biopython_to_mmcif import BiopythonToMmcifResidueIds, ResidueId, BioResidueId


class SSType(Enum):
    HELIX = auto()
    STRAND = auto()
    NO_TYPE = auto()  # no ss type found for a residue


def contained_in_segment(number: int, segment_starts: List[int], segment_ends: List[int]):
    """ Returns True if number is contained in any of the non-overlapping segments. Lists have to be sorted. """
    i = -1 + bisect(segment_starts, number)  # bisect: a[:i] have e <= x, and all e in a[i:] have e > x
    return i >= 0 and number <= segment_ends[i]


@dataclass
class SSForChain:
    """ Can determine if a residue has a certain secondary structure by bisecting the ss_start lists.

    <ss>_start and <ss>_end lists have to be sorted and for a ss type their positions correspond (e.g. helices_start[i] and helices_end[i]
    mark the start and end of i-th helical segment)

    If we were always to compare the whole sequences of residues, it would be faster to implement it with a pointer to
    current SS start/end segment. But this should be fast enough, is easier to implement (just using bisect module),
    and more versatile use (fast for single residues).
    """
    helices_start: List[int]
    helices_end: List[int]
    strands_start: List[int]
    strands_end: List[int]

    def _is_ss_for_residue(self, residue_serial: int, ss_start: List[int], ss_end: List[int]) -> bool:
        return contained_in_segment(residue_serial, ss_start, ss_end)

    def ss_for_residue(self, r: ResidueId):
        residue_serial = r.label_seq_id  #  davam author_seq_id nekam, kde mam mit label_seq (nebo to predelej na author, ale vykaslat se na insertion code)

        for ss_type, ss_start, ss_end in ((SSType.HELIX, self.helices_start, self.helices_end), (SSType.STRAND, self.strands_start, self.strands_end)):
            if self._is_ss_for_residue(residue_serial, ss_start, ss_end):
                return ss_type

        return SSType.NO_TYPE


@dataclass
class SSForStructure:
    ss_for_chains: Dict[str, SSForChain]

    def ss_for_residue(self, r: ResidueId):
        return self.ss_for_chains[r.chain_id].ss_for_residue(r)


TResidueData = TypeVar('TResidueData')


@dataclass(frozen=True, eq=False)
class SetOfResidueData(Generic[TResidueData]):
    # TResidueId = Tuple[int, str]  # author_seq_id, author_insertion_code? As in BioPython, but no hetero flag

    data: List[TResidueData]
    structure_id: str

    def get_full_id(self):
        raise NotImplementedError

    def __len__(self):
        return len(self.data)

    def __iter__(self) -> Iterator[TResidueData]:
        return iter(self.data)

    def __getitem__(self, item):
        return self.data[item]

    def __bool__(self):
        return bool(self.data)

    def serialize(self):
        return self.get_full_id()

    def __add__(self, other: 'SetOfResidueData[TResidueData]'):
        return CombinedSetOfResidueData[TResidueData].create_from_groups_of_residues(self, other)

    def __hash__(self):
        return hash(self.get_full_id())  # zatim asi hack, do budoucna trochu jinak bude cachovani fungovat (tohle by stejně asi nefunguvalo -- kdyby tam sla instance, nemohl by to GC odstranit

    def __eq__(self, other):
        return self.get_full_id() == other.get_full_id()

    @staticmethod
    def residues_from_label_seq_ids(label_seq_ids: Iterable[int], mapping: BiopythonToMmcifResidueIds.Mapping, bio_chain: Chain):
        return list(
            map(lambda bio_res_id: bio_chain[bio_res_id],
                # get a list of Biopython Residue ids
                map(lambda label_seq_id: mapping.to_bio(label_seq_id), label_seq_ids)
            )
        )

    def get_atoms(self: 'SetOfResidueData[Residue]', exclude_H=True) -> Iterable[Atom]:
        """ exclude_H by default, so the method is compatible with freesasa, see more in `get_sasa` fn

        in Python one cannot specialize SetOfResidueData for SetOfResidueData[Residue] and provide its own definition,
        like in c++, so this is the way how to do it. Note the type annotation at `self` argument."""

        for r in self:
            for a in r.get_atoms():
                # optionally exclude hydrogen atoms
                if exclude_H and a.element == 'H':
                    continue

                yield a


@dataclass(frozen=True, eq=False)
class CombinedSetOfResidueData(SetOfResidueData[TResidueData]):
    combined_id: Tuple[Any, ...]

    @classmethod
    def create_from_groups_of_residues(cls, *groups_of_residues: SetOfResidueData[TResidueData]):
        all_residues = []

        for residues in groups_of_residues:
            all_residues.extend(residues)

        return CombinedSetOfResidueData[TResidueData](all_residues,
                                     groups_of_residues[0].structure_id,
                                     tuple(sorted(g.get_full_id() for g in groups_of_residues)))  # sorted, so __eq__ works as expected

    def get_full_id(self):
        return self.combined_id


@dataclass(frozen=True, eq=False)
class ChainResidueData(SetOfResidueData[TResidueData]):
    chain_id: str

    def get_full_id(self):
        return (self.structure_id, self.chain_id)


@dataclass(frozen=True, eq=False)
class DomainResidueData(SetOfResidueData[TResidueData]):
    chain_id: str
    domain_id: str

    def get_full_id(self):
        return (self.structure_id, self.chain_id, self.domain_id)


class SetOfResidues(SetOfResidueData[Residue]):
    pass


class ChainResidues(ChainResidueData[Residue], SetOfResidues):
    @classmethod
    def from_chain(cls, chain: Chain, mapping: BiopythonToMmcifResidueIds.Mapping):
        return cls(
           [chain[bio_id] for bio_id in mapping.label_seq_id__to__bio_pdb.values()],
            chain.get_parent().get_parent().id,
            chain.id,
        )


@dataclass
class DomainResidueMapping:
    domain_id: str
    chain_id: str
    segment_beginnings: List[int]
    segment_ends: List[int]

    @classmethod
    def from_domain_on_another_chain(cls, domain: 'DomainResidueMapping', new_domain_chain_id: str,
                                     label_seq_id_offset: int):

        # todo nebo bych mohl vybrat třeba pomocí těch segmentů (a +1 garance by nebyla třeba), kdybych mohl
        #  indexovat od do label seq id, ale to by bylo ordered dict, jestli to vubec jde.
        return cls(domain.domain_id, new_domain_chain_id,
                   [i + label_seq_id_offset for i in domain.segment_beginnings],
                   [i + label_seq_id_offset for i in domain.segment_ends],
       )

    def __iter__(self) -> Iterator[int]:
        """ Returns a sequence of residue label_seq_id in a domain """
        # todo again here I suppose label_seq_id is +1 sequential..
        for segment_start, segment_end in zip(self.segment_beginnings, self.segment_ends):
            yield from range(segment_start, segment_end + 1)

    def __contains__(self, label_seq_id):
        return contained_in_segment(label_seq_id, self.segment_beginnings, self.segment_ends)

    def __len__(self):
        residue_count = 0
        for segment_start, segment_end in zip(self.segment_beginnings, self.segment_ends):
            residue_count += segment_end + 1 - segment_start

        return residue_count

    def to_set_of_residue_ids(self, structure_id: str, skip_auth_seq_id=lambda id: False) -> SetOfResidueData[ResidueId]:
        return DomainResidueData[ResidueId]([ResidueId(auth_seq_id, ' ', self.chain_id) for auth_seq_id in self
                                             if not skip_auth_seq_id(auth_seq_id)], structure_id, self.chain_id, self.domain_id)

    def get_spans(self):
        return list(zip(self.segment_beginnings, self.segment_ends))


# todo vlastne byl problem, ze jsem nekdy dostal domén víc, takže jsem je nevěděl, jak namapovat mezi sebou
# uvažoval jsem nějaký podrozdělování domén na víc a tak podobně, to je blbost
#  - přeskupovat segmenty do nových domén je blbost
# Myslim, že stačí pro každou doménu najít tu s největším overlappem (největší součet velikostí merged_segmentů)
# což je n^2, ale kolik je max domén? To zatim zjistit nemůžu, ale asi moc ne, když max len lcs je ...
# navíc vyberu jen takový SEGMENTY, co maj překryv s lcs, abych to urychlil?
# max lcs len je 4620 a median 127
# jako mohl bych ty domeny nejak seradit, ale oni obecne muzou skakat přes sebe, tak nevim, jestli je něco rychlejšího
# než D^2?
# v čase n*D můžeš říct na kterých pozicích které domény
# domény jsou disjunktní (hopefully), to by mohlo pomoct ve složitosti
#   tzn. ke každý pozici budu mít max 1 apo a 1 holo
#   - to mi dá vlastně possible páry
#   - bral bych jaccard similarity (IoU)
#   - jasne no, párů může být O(min(D^2, N))
#   jejich ohodnocení zvládnu naivním algoritmem v čase Theta(D^2 * mean_seg_count)
#    - a lepším algoritmem O(mean_seg_count * D), kdy spočítám velikosti domén + při passu přes N (může být segment-level, nikoliv res-level) spočítám
#    - ten lepší by měl list segmentů v apo, pak holo plus domén původu pro každý ten segment
#           - pak by mergoval segmenty úplně stejně, jako v merge_segments (+ poznamenávat, z kt. 2 domén jsou)
#           - pak projít ty merged_segmenty a do dictu (d1, d2) -> int přičítat vždy velikost segmentu těch dvou domén
#           - nakonec pro všechny keys (d1, d2) toho dictu spočítat jaccard? similarity.

def merge_segments(segments1, segments2):
    """ For two arrays of segments, do an AND-like operation. That is return segments that describe indices that are
    in both input segments.

    :param segments1: array-like; list of tuples (start/end) or np.ndarray of shape (n_segments, 2)
    :param segments2: same as above
    :return: np.ndarray with shape (n_result se
    """
    segments1, segments2 = map(np.asarray, (segments1, segments2))
    common_segments = []

    i = j = 0
    while i < len(segments1) and j < len(segments2):
        s1_beg, s2_beg = segments1[i, 0], segments2[j, 0]
        s1_end, s2_end = segments1[i, 1], segments2[j, 1]

        maybe_seg_beg = max(s1_beg, s2_beg)
        maybe_seg_end = min(s1_end, s2_end)

        if maybe_seg_beg <= maybe_seg_end:
            common_segments.append((maybe_seg_beg, maybe_seg_end))

        # if a segment was fully consumed, move to next one
        if s1_end == maybe_seg_end:
            i += 1
        if s2_end == maybe_seg_end:
            j += 1

    return np.array(common_segments)


def merge_domains(d1: DomainResidueMapping, d2: DomainResidueMapping, label_seq_offset: int):
    d1_segments = np.stack([d1.segment_beginnings, d1.segment_ends], axis=-1)
    d2_segments = np.stack([d2.segment_beginnings, d2.segment_ends], axis=-1) - label_seq_offset

    common_segments = merge_segments(d1_segments, d2_segments)

    common_segment_begs, common_segment_ends = np.array(common_segments).T

    new_d1 = DomainResidueMapping(d1.domain_id, d1.chain_id,
                                  common_segment_begs.tolist(),
                                  common_segment_ends.tolist())

    new_d2 = DomainResidueMapping(d2.domain_id, d2.chain_id,
                                  (common_segment_begs + label_seq_offset).tolist(),
                                  (common_segment_ends + label_seq_offset).tolist())
    return new_d1, new_d2


class DomainResidues(DomainResidueData[Residue], SetOfResidues):
    # todo pošlu sem chainresidues a jejich label_seq_id (chci jen observed v obou sekvencich, podvybrat domenu)
    @classmethod
    def from_domain(cls, domain: DomainResidueMapping, bio_structure: Model,
                    residue_id_mapping: BiopythonToMmcifResidueIds.Mapping, skip_label_seq_id=lambda id: False):
        bio_chain = bio_structure[domain.chain_id]  # todo proč??? proč tam nepošlu rovnou chain?

        domain_residues = [bio_chain[residue_id_mapping.to_bio(label_seq_id)]
                           for label_seq_id in domain if not skip_label_seq_id(label_seq_id)]

        return cls(domain_residues, bio_structure.get_parent().id, domain.chain_id, domain.domain_id)

    def get_spans(self, residue_id_mapping: BiopythonToMmcifResidueIds.Mapping, auth_seq_id=False):
        """Get domain spans (segments).

        :param residue_id_mapping:
        :param auth_seq_id: if True, return auth_seq_ids instead of label_seq_ids (however, it works just with integer
        ids, ignores insertion codes)
        :return: Domain segments. List of tuples (from, to_inclusive), of label_seq_ids (or auth_seq_ids)
        """

        # data is sorted (because DomainResidueMapping is), so just break a span after label_seq_id not contiguous
        def get_res_id(r):
            if auth_seq_id:
                return r.id[1]
            return residue_id_mapping.to_label_seq_id(BioResidueId(*r.id))

        try:
            residues = iter(self.data)

            span_start = get_res_id(next(residues))
            last_label_seq_id = span_start

        except StopIteration:
            return []

        spans = []

        for r in residues:
            id_ = get_res_id(r)
            if id_ - last_label_seq_id > 1:
                # new span, if label_seq_ids weren't contiguous
                spans.append((span_start, last_label_seq_id))
                span_start = id_

            last_label_seq_id = id_

        spans.append((span_start, last_label_seq_id))

        return spans


@dataclass
class ScrewMotionResult:
    """ Represents result of GetHingeAngle

    Angle in radians, translations in angstroms, all absolute values."""
    angle: float
    translation_in_axis: float
    translation_overall: float
