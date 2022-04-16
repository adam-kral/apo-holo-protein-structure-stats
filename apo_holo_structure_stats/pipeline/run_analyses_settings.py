import functools
import multiprocessing
import shelve
from functools import partial, lru_cache
from pathlib import Path

from apo_holo_structure_stats.core.analyses import *
from apo_holo_structure_stats.core.base_analyses import Analyzer
from apo_holo_structure_stats.core.mp_lru_cache import MultiProcessingLRUCache as mulproc_lru_cache


# todo in future the dependecy instantiation might be automatic, one would just initialize the variables with
#  lambdas that would wrap the created instances?

# todo wrap in simple memoize cache or API requests in MultiProcessingLRUCache
# lightweight cache is only in Python 3.9, for 3.8 need to use lru_cache even for maxsize=None
from apo_holo_structure_stats.input.download import parse_mmcif


def simplify_arg(arg):
    try:
        # works with SetOfResidueData, so that in memory there can be only its id left,
        # not the whole BioPython's structure possibly..
        return arg.get_full_id()
    except AttributeError:
        return arg

mulproc_lru_cache = partial(mulproc_lru_cache, args_to_key_fn=simplify_arg)


def simplify_args_decorator(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        args = map(simplify_arg, args)
        kwargs = {k: simplify_arg(v) for k, v in kwargs.items()}
        return func(*args, **kwargs)
    return wrapper


# nejde writovat ve více procesech/vláknech
# read je ale ok
# takze read je jednoduchy, misto get_ss to definuju jako load result (asi teda ta shelf picklovana? Nebo gdbm se stringem)

# ukladani udelam pres iterovani futurů? Klidne fast open a once in a while zavolam sync..

# jediny, co si musim pamatovat je name db - to je normalni, input output...
# jak bude velika zhruba ta db? Asi ji nebudu shardovat do jobu, proste ji tam hodim snad celou (treba vetsi joby udelam...)

# takze to bude na predvypocet api veci, s danym limitem? A co kdybych chtěl i nacitat ty struktury?
# nebo jiny (rychlejsi api bez omezeni)
# to ted resit nebudu, kazdej si muze udelat vlastni meziskript...
# ale bylo by hezky, kdyby slo pridavat single struct level analyzy nejak.. (ne-api)
# jako treba by to slo do filter structures (a dat tam i to getsasa)

# chci ty single struct veci (krome nutneho is_holo) delat jenom pro proteiny v pairs,
# podle me to budou ale skoro vsechny stejne...


class SavedAnalysis(Analyzer):
    def __init__(self, db_filename, analysis):
        self.db_filename = db_filename
        self.db = None

        self.name = analysis.get_name()

        super().__init__()

    def run(self, key):
        # not opening db file in init, because it wouldn't be picklable (to other processes)
        if self.db is None:
            self.db = shelve.open(self.db_filename, flag='r')

        try:
            return self.db[key]
        except KeyError:
            raise MissingDataException

    def get_name(self):
        return self.name






# todo kolik je ne páru v jedný skupině, ale kolik je tříd ekvivalence/tranzitivních párů
MP_CACHE_SIZE = 1e4

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

"""
vsechny cached veci jsou blbost celkem.. Chtel bych si je ukladat na disk tak jako tak (domains treba, ne nutne 
ss nebo sasa). Navic ze sasa nema cenu cachovat to s1+s2 preci... To ale specifikovat uplne nejde.. Treba arg no_cache
v GetInterfaceBuriedArea.run kdyz volam get_sasa? Moc clean to neni, navic by to spadlo, pokud bych to neobalil v lru_cache
ktera by ten arg vyzobla..

Navic bych nepretizil servery predvypoctem (i kdyz sasa - tam se hodi asi to cachovani, protoze uz musim nacist ty struktury)

"""


# todo on end do SavedAnalysis.db.close


# todo delete this
#  need to skip non-xray now (we still have old data unskipped in `filter_structures`)
class NotXrayDiffraction(Exception):
    pass

def parse_mmcif_exp_method_hack(not_xray: dict, pdb_code: str):
    ''' as we want to cache the parsed structures (we do 2), we also want to record that the structure is not x-ray,
    so we need to parse it only once. (this is a hack see above, is solved in `filter_structures` '''
    # not_xray a dict, as multiprocessing.Manager has not .set() attr
    if pdb_code in not_xray:
        raise NotXrayDiffraction

    apo_parsed, apo_metadata = parse_mmcif(pdb_code, with_extra=True, allow_download=False)
    # apo_parsed, apo_metadata = parse_mmcif(pdb_code, with_extra=True, allow_download=True)  # todo temp hack for testing

    # todo delete this
    #  need to skip non-xray now (we still have old data unskipped in `filter_structures`)
    if apo_metadata.mmcif_dict['_exptl.method'][0] != 'X-RAY DIFFRACTION':
        not_xray[pdb_code] = True
        raise NotXrayDiffraction(f"method is {apo_metadata.mmcif_dict['_exptl.method'][0]}")

    return apo_parsed


def configure_pipeline(manager, input_dir: Path):
    # todo all other analyses could be cached at least with lru_cache for the time a pair is processed (they are called
    # multiple times in some cases..)
    # def end_pair_callback():
    #     for cache in caches:
    #         cache.cache_clear()


    # todo delete this
    #  need to skip non-xray now (we still have old data unskipped in `filter_structures`)
    not_xray = manager.dict()
    parse_mmcif = mulproc_lru_cache(partial(parse_mmcif_exp_method_hack, not_xray), manager, max_size=2)
    # parse_mmcif = mulproc_lru_cache(partial(parse_mmcif, allow_download=False), manager, max_size=2)


    get_chains = GetChains()

    get_c_alpha_coords = GetCAlphaCoords()
    get_centroid = GetCentroid((get_c_alpha_coords,))
    get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
    get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))

    get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))
    get_rmsd = GetRMSD((get_centered_c_alpha_coords, get_rotation_matrix))

    # get_ss = GetSecondaryStructureForStructure()
    # get_ss = mulproc_lru_cache(get_ss, manager, MP_CACHE_SIZE)
    get_ss = SavedAnalysis(str(input_dir / 'db_get_ss'), GetSecondaryStructureForStructure())
    cmp_ss = CompareSecondaryStructure((get_ss,))

    get_domains = SavedAnalysis(str(input_dir / 'db_get_domains'), GetDomainsForStructure())

    get_sasa = GetSASAForStructure()
    # nepouzivat cache, protoze hranice domen se mohou menit podle partnera v páru (unobserved nebo cropped protoze
    # lcs)
    # get_sasa = mulproc_lru_cache(get_sasa, manager, MP_CACHE_SIZE)
    get_interdomain_surface = GetInterfaceBuriedArea((get_sasa,))

    comparators__residues_param = [get_rmsd]
    comparators__residue_ids_param = [cmp_ss]

    comparators__domains__residues_param = [get_rmsd]
    comparators__domains__residue_ids_param = [cmp_ss]

    comparators__2DA__residues_param = [get_hinge_angle]

    analyzer_namespace = locals()
    del analyzer_namespace['manager']  # remove the argument from the namespace
    del analyzer_namespace['input_dir']  # remove the argument from the namespace

    # return dotdict(analyzer_namespace)  # najednou chyba, kdyz jsem to spoustel vedle v cryptic_binding_sites, tak to
    # radsi opravim i tu (navic to nakonec stejne bylo zbytecny)
    return analyzer_namespace
