from functools import lru_cache
from typing import NamedTuple, Dict, Tuple


class Polymer:
    pass

# domain, chain, subsets of polymers

class Ligand:
    pass
# get_atoms


# class AnalysisArgs:
#     # bud RAII nebo metoda run
#     def __init__(self, *args, **kwargs):
#         pass
#
#     def serialize(self):
#         """ Returns analysis name and serialized arguments """
#         pass
#
#
# class AnalysisResult:
#     # tohle spis bude metoda na Analyzer -> serializeResult
#     # protoze nechci vzdycky delat Analyzer()(args).value ...
#
#     def serialize(self):
#         pass
#     #
#     # def value(self):
#     #     pass


class Analyzer:
    """ Instance behaves like a function.

    Use analyzer_instance(arg1, arg2) to get the result. """
    dependencies: Tuple['Analyzer', ...]

    def __init__(self, dependencies: Tuple['Analyzer', ...] = ()):
        self.dependencies = dependencies

    def __call__(self, *args, **kwargs):
        """ Use this """
        return self.run(*args, *self.dependencies, **kwargs)

    def run(self, *args, **kwargs):
        """ Provide your own method. """
        raise NotImplementedError

    def get_name(self):
        return type(self).__name__


class CachedAnalyzer(Analyzer):
    """ Holds instances of analyses and their results. Analyses are run though this object. But then -- parallelization? Serialization?"""

    def __call__(self, *args, **kwargs):
        return self.analysis_cached_run(*args, **kwargs)

    # @lru_cache(maxsize=None)
    def analysis_cached_run(self, *args, **kwargs):
        return super().__call__(*args, **kwargs)


class SerializableAnalyzer(Analyzer):
    def serialize(self, result, *args, **kwargs):
        """ default implementation assumes result, args and kwargs are primitive types, therefore know how to serialize to a string value.

        :param result: result of __call__/run
        :param args: for now: arbitrary args, serve as identifier for an analysis result in the serialized data. (Don't need to be actual args to __call__/run)
        :param kwargs: similar as args
        :return: a potentially nested dict/list of primitive types
        """

        def serialize_argument(a):
            try:
                return a.serialize()  # třeba takhle přidat custom argumenty, ktery se budou umet serializovat? Mozna to jeste muzu ujasnit pomoci typing?
            except AttributeError:
                return a

        serialized_args = [serialize_argument(a) for a in args]
        serialized_kwargs = {k: serialize_argument(v) for k, v in kwargs.items()}

        return {'analysis_name': self.get_name(), 'args': serialized_args, 'kwargs': serialized_kwargs, 'result': result}  # jeste serialize result?



class SerializableCachedAnalyzer(CachedAnalyzer, SerializableAnalyzer):
    pass



# do budoucna -- mohu udělat dekorátor, který z funkce udělá  Analyzer classu, (dobrý nápad??), dependencies a base classa jako parametry decorátoru, funkce pak bude brát ty dependecies
# spíš než dekorátor možná

# serializace - non-output mohou mít i třeba pickle?
#      - volá se při prvním výsledku analýzy (v třídě Analyzer). Serializace teoreticky může být do různých formátů, mám ale json (a možná pickle)
# jednoznačná identifikace argumentů??
    # is - equality? (nebude fungovat po un-picklování file-cached dat)
# Jednotlivé analyzery -- depedency injection?
# definovat - instance do proměnných a proměnné do dependencies i celkového seznamu analyzátorů (projde se a vykoná analýzy)
# zatim cached podle __hash__ argumentů, tzn. povětšinou in-memory identity objektu (do budoucna lze předefinovat __hash__ na contents AnalysisResult
#
#
# @wrap_in_analyzer(CachedAnalyzer)
# def get_ligands(struct):
#     return itertools.chain(get_hetero_atom_residues(struct), get_short_peptide_ligands(struct, 15))
#
