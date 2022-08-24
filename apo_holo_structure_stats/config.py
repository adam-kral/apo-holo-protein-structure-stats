""" Not used, but one of the other (simpler) but worse options (no IDE autocomplete/mypy type check). """
import os
import logging
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)


class ConfigNotFoundException(Exception):
    pass


if os.environ.get('AH_CONFIG'):
    try:
        config = yaml.safe_load(Path(os.environ['AH_CONFIG']))
    except FileNotFoundError as e:
        raise ConfigNotFoundException(f'Config `{os.environ["AH_CONFIG"]}` specified in AH_CONFIG '
                                      f'environment variable not found') from e
else:
    config = yaml.safe_load(Path(__file__).parent / 'sample_config.yaml')
    logger.warning(f'Config not specified in AH_CONFIG environment variable. Using `sample_config.yaml`.')
