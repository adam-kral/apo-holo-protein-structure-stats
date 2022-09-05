import logging
from argparse import ArgumentParser
from os import environ
from pathlib import Path

import yaml


# def get_argument_parser(add_config=True):
def get_argument_parser(**kwargs):
    parser = ArgumentParser(**kwargs)
    add_loglevel_args(parser)
    #
    # if add_config:
    #     parser.add_argument('-c', '--config', type=Path)

    return parser

def add_loglevel_args(parser: ArgumentParser):
    parser.add_argument(
        '--debug',
        action="store_const", dest="loglevel", const=logging.DEBUG,
        default=logging.WARNING,
    )
    parser.add_argument(
        '-v', '--verbose',
        action="store_const", dest="loglevel", const=logging.INFO,
    )


def get_config(path: Path = None):
    default_config = Path(__file__, '..', '..', '..') / 'default_config.yaml'
    config = Path(environ.get('AH_CONFIG', default_config))

    if path:
        config = path

    return yaml.safe_load(config)
