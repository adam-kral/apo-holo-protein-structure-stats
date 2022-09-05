import functools
import multiprocessing
import shelve
from functools import partial, lru_cache
from pathlib import Path

from apo_holo_structure_stats.core.analyses import *
from apo_holo_structure_stats.core.base_analyses import Analyzer
from apo_holo_structure_stats.core.mp_lru_cache import MultiProcessingLRUCache as mulproc_lru_cache

from apo_holo_structure_stats.input.download import parse_mmcif as parse_mmcif_fn


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
# je jich o dost min, tim spis, pokud nas bude zajimat nejaky uzsi dataset nez uplne vsechny struktury z PDB


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

"""
(skoro) vsechny cached veci jsou blbost celkem.. Chtel bych si je ukladat na disk tak jako tak (domains treba, ne nutne 
ss nebo sasa). Navic ze sasa nema cenu cachovat to s1+s2 preci... To ale specifikovat uplne nejde.. Treba arg no_cache
v GetInterfaceBuriedArea.run kdyz volam get_sasa? Moc clean to neni, navic by to spadlo, pokud bych to neobalil v lru_cache
ktera by ten arg vyzobla..

"""


# todo on end do SavedAnalysis.db.close


def configure_pipeline(input_dir: Path, manager=None):
    # todo all other analyses could be cached at least with lru_cache for the time a pair is processed (they are called
    # multiple times in some cases..)
    # def end_pair_callback():
    #     for cache in caches:
    #         cache.cache_clear()

    # todo multiprocessing cache is, especially for large objects such as Bio.Structures, counter-productive due to slow IPC,
    # (e.g. had 4 processors, but only 1.5 were used on avg).
    # do not use multiprocessing, or use better strategy - larger chunks than one pair,
    # parse_mmcif = mulproc_lru_cache(oarse_mmcif_fn, manager, max_size=2)

    parse_mmcif = lru_cache(maxsize=3)(parse_mmcif_fn)  # this cache is important, saves most of the  re-loading (large)
    # structures, especially useful as large structures have many chains (and pairs) and would be re-loaded
    # disproportionately often

    get_chains = GetChains()

    get_c_alpha_coords = GetCAlphaCoords()
    get_centroid = GetCentroid((get_c_alpha_coords,))
    get_centered_c_alpha_coords = GetCenteredCAlphaCoords((get_c_alpha_coords, get_centroid))
    get_rotation_matrix = GetRotationMatrix((get_centered_c_alpha_coords,))

    get_hinge_angle = GetHingeAngle((get_c_alpha_coords, get_centroid, get_rotation_matrix))
    get_rmsd = GetRMSD((get_centered_c_alpha_coords, get_rotation_matrix))

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

    return analyzer_namespace
