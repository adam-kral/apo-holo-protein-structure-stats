import logging
from argparse import ArgumentParser


def add_loglevel_args(parser: ArgumentParser):
    parser.add_argument(
        '-d', '--debug',
        action="store_const", dest="loglevel", const=logging.DEBUG,
        default=logging.WARNING,
    )
    parser.add_argument(
        '-v', '--verbose',
        action="store_const", dest="loglevel", const=logging.INFO,
    )
