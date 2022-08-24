""" Defines default settings and imports possibly user-defined local setting, overwriting some of the defaults

Local settings can overwrite the values of these settings like this:
    KEY = NEW_VAL
or like this:
    from apo_holo_structure_stats.settings import klass
    klass.KEY = NEW_VAL

Path user-defined local setting if supplied via the enviroment

Local settings can also contain additional settings for user-created code. In that case you can import them through
`apo_holo_structure_stats.settings` or perhaps more conveniently directly from your local settings module (you'll have
the IDE code completion).
"""
import importlib
import sys
from pathlib import Path
import os
from typing import Callable

import yaml

STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path(os.environ.get('STRUCTURE_DOWNLOAD_ROOT_DIRECTORY', 'pdb_structs'))

# todo for testing, probably create an env var.. IS_PRODUCTION=1 or something like that
# STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path(Path(__file__).parent, 'paper_repl', 'pdb_structs')

MIN_STRUCTURE_RESOLUTION = 2.5

class LigandSpec:
    MIN_NON_H_ATOMS = 6
    MIN_RESIDUES_WITHIN_LIGAND = 6
    MIN_RESIDUES_WITHIN_LIGAND__RADIUS = 4.5


import yaml

# @dataclasses.dataclass
# class BaseConfig(yaml.YAMLObject):
#     name: str
#     items: List[str] = field(default_factory=lambda: ['a', 'b', 'c'])
#     yaml_tag: str = u'!Yaz'


def import_local_settings_module_from_file(path: str):
    import importlib.util
    import sys

    module_name = 'local_settings'  # unimportant, just placeholder, but hopefully doesn't need to match file name

    # https://docs.python.org/3/library/importlib.html#importing-a-source-file-directly
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)


class LocalSettingsNotFoundException(Exception):
    pass


# possibly import local settings
file_path_or_module_path = os.environ.get('AH_LOCAL_SETTINGS')

if file_path_or_module_path:
    if file_path_or_module_path.endswith('.py'):
        # path/to/my/module.py
        try:
            import_local_settings_module_from_file(file_path_or_module_path)
        except FileNotFoundError as e:
            raise LocalSettingsNotFoundException(f'File `{file_path_or_module_path}` not found. '
                                                 f'\n\nIf you are using also _your code_ and not just modifying the settings'
                                                 f' you can specify just the module name such as `my_module` or '
                                                 f'`package.subpackage.my_module` - maybe you accidentally '
                                                 f'appended .py to the module name?') from e

    else:
        # my_package.my_module
        try:
            importlib.import_module(file_path_or_module_path)
        except ImportError as e:
            raise LocalSettingsNotFoundException(f'Module `{file_path_or_module_path}` not found. Provide module name, '
                                                 f'possibly including the packages it is located in.') from e


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
    """
    STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path('pdb_structs')
    MIN_STRUCTURE_RESOLUTION = 2.5

    class LigandSpec:
        MIN_NON_H_ATOMS = 6
        MIN_RESIDUES_WITHIN_LIGAND = 6
        MIN_RESIDUES_WITHIN_LIGAND__RADIUS = 4.5

    API_REQUESTS_TIMEOUT = 10
    API_REQUESTS_RETRIES = 5


def from_nested_dict(klass: type, d: dict) -> None:
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
                raise Exception from e

            used_keys.add(key)

    if used_keys != d.keys():
        raise Exception


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

#
# print(to_nested_dict(Settings))
#
print(yaml.safe_dump(to_nested_dict(Settings)))
#
# # from_nested_dict(Settings,
# # {'yaml_tag': 'hovno', 'LigandSpec': {'yaml_tag': 'sdfa'}}
# # )


yaml_settings = os.environ.get('AH_SETTINGS_FILE')
if yaml_settings:
    # get the settings class (default is the above Settings, but can be its descendant,
    #   should be directly in a module and not inside any other class namespace)
    parts = os.environ.get('AH_SETTINGS_CLASS', 'apo_holo_structure_stats.settings.Settings').rsplit('.', 1)
    if len(parts) != 2:
        sys.stderr.write(f'AH_SETTINGS_CLASS should be in format `[packages.]module.Class`. Received `{".".join(parts)}`')
        sys.exit(1)

    module_name, class_name = parts
    module = importlib.import_module(module_name)
    settings_class = getattr(module, class_name)

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
