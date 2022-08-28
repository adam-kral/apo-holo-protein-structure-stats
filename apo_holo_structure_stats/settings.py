""" Defines the default settings and possibly overwrites them using local settings from yaml (specified
in AH_SETTINGS_FILE environment variable). Additionally, some settings are loaded directly from envvars.

Individual settings are in the Settings class. (That has some disadvantages - if it were in a module,
classes inside the module could be imported, as opposed to classes nested in a class. But also some advantages
- easier implementation of serialization to yaml, for extended settings.)
This handling of Settings allows code completion and mypy type checking.

Running this module directly (as __main__) prints (the attributes of) Settings class serialized in yaml.
(Or that of AH_SETTINGS_CLASS, see Settings' docstring.) Useful in generating the default config.

The yaml is required to specify _all_ the attributes of the Settings class.
"""
import importlib
import sys
from pathlib import Path
import os
from typing import Callable

import yaml

# 1) define Settings with defaults


class Settings:
    """ Settings of the framework. Can be loaded from/dumped to yaml.

    Specify the yaml file path in the environment variable AH_SETTINGS_FILE.

    If you want to add more settings, and want them to be in the same single yaml, you can inherit this class
    and add more class attributes. Then set envvar AH_SETTINGS_CLASS.

    Classes can be nested, only supported attribute types are those supported by pyyaml and pathlib.Path.

    Attributes starting with `_` are skipped and all callables (except classes = namespaces) too.

    (De)serialization is handled by `from_nested_dict` and `to_nested_dict` functions below.

    The yaml has to contain all the (serialized) attributes of the class (including its superclasses). After
    loading the yaml, everything will be updated, incl. its superclasses.

    This class is not meant to be instantiated.
    """
    STRUCTURE_STORAGE_DIRECTORY = Path('pdb_structs')
    MIN_STRUCTURE_RESOLUTION = 2.5

    class LigandSpec:
        MIN_NON_H_ATOMS = 6
        MIN_RESIDUES_WITHIN_LIGAND = 6
        MIN_RESIDUES_WITHIN_LIGAND__RADIUS = 4.5
        PEPTIDE_MAX_LENGTH = 15

    MIN_OBSERVED_RESIDUES_FOR_CHAIN = 50

    API_REQUESTS_TIMEOUT = 10
    API_REQUESTS_RETRIES = 5

    FILTER_STRUCTURES_CLASS = 'apo_holo_structure_stats.pipeline.filter_structures.StructureProcessor'
    MAKE_PAIRS_CLASS = 'apo_holo_structure_stats.pipeline.make_pairs_lcs.ApoHoloMatchmaker'


# 2) load settings from yaml (possibly)

class ClassSerializerException(Exception):
    pass


def from_nested_dict(klass: type, d: dict) -> None:
    """ Converts class to nested dict (attr_name -> attr_value). The classes can be nested.

    Skips attributes starting with `_`, incl. magic attrs. Skips callables, except the classes.
    All attributes need to be provided in the dict. """
    used_keys = set()

    for c in klass.__mro__:
        for key, value in c.__dict__.items():
            # ignore magic methods or anything starting with '_'
            if key.startswith('_'):  # or I could consider only UPPERCASE, just like in django settings
                continue

            if isinstance(value, Callable) and not isinstance(value, type):
                # skip functions/methods (but not classes = useful namespaces)
                continue

            try:
                if isinstance(value, type):
                    # attribute is a class -> load it from the nested dict recursively
                    from_nested_dict(value, d[key])
                elif isinstance(value, Path):
                    setattr(c, key, Path(d[key]))  # create Path from string d[key]
                else:
                    setattr(c, key, d[key])
            except KeyError as e:
                raise ClassSerializerException(f'Missing key (class attribute) {key} in the argument `d` for class '
                                               f'{klass}') from e

            used_keys.add(key)

    if used_keys != d.keys():
        raise ClassSerializerException(f'Too many keys provided in the argument `d` for class {klass}.'
                                       f'\n Used keys: {used_keys}. Provided keys: {set(d.keys())}')


def to_nested_dict(klass: type) -> dict:
    d = {}
    for c in reversed(klass.__mro__):  # reversed so that attrs that are shadowed in superclasses are correctly overwritten
        for key, value in c.__dict__.items():
            if key.startswith('_'):
                continue

            if isinstance(value, Callable) and not isinstance(value, type):
                # skip functions/methods (but not classes = useful namespaces)
                continue

            if isinstance(value, type):
                # attribute is a class -> convert it to the nested dict recursively
                value = to_nested_dict(value)
            elif isinstance(value, Path):
                value = str(value)
            d[key] = value

    return d


def load_class_from_str(class_location: str):
    parts = class_location.rsplit('.', 1)
    if len(parts) != 2:
        raise ValueError('class_location should be in format `[packages.]module.Class')

    module_name, class_name = parts
    module = importlib.import_module(module_name)
    return getattr(module, class_name)


# run this module directly to print the default settings
if __name__ == '__main__':
    print(yaml.safe_dump(to_nested_dict(Settings), sort_keys=False))


# get the settings class (default is the above Settings, but can be its descendant,
#   should be directly in a module and not inside any other class namespace)
settings_class_location = os.environ.get('AH_SETTINGS_CLASS', f'{__name__}.Settings')
try:
    settings_class = load_class_from_str(settings_class_location)
except ValueError:
    sys.stderr.write(f'AH_SETTINGS_CLASS should be in format `[packages.]module.Class`. Received `{settings_class_location}`')
    sys.exit(1)

# load run-specific settings, if supplied
yaml_settings = os.environ.get('AH_SETTINGS_FILE')
if yaml_settings:
    # load settings class attributes from the yaml settings dictionary (will replace the default ones in the code)
    try:
        with open(yaml_settings) as f:
            settings_dict = yaml.safe_load(f)
    except FileNotFoundError as e:
        sys.stderr.write(f'You specified yaml settings in `AH_SETTINGS_FILE` environment variable '
                         f'but the settings file `{yaml_settings}` was not found.')
        sys.exit(1)
    print(settings_dict)
    from_nested_dict(settings_class, settings_dict)

# 3) some settings are loaded directly from envvars
# replace `STRUCTURE_STORAGE_DIRECTORY` setting if explicitly in envvar
if dir_str := os.environ.get('AH_STRUCTURE_STORAGE_DIRECTORY'):
    Settings.STRUCTURE_STORAGE_DIRECTORY = Path(dir_str)