"""
Automatically runs the pipeline scripts. (If a test case class is run, it will run all of the pipeline scripts in serial on
two datasets (subsampled pdb) or the paper dataset.) Nothing is mocked, i.e. it downloads structures from pdb and makes requests
to APIs. This is meant to be for local testing right before running it in batches, e.g. on Metacentrum where debugging is not
desirable.

The output directories are not deleted (so that you can rerun specific tests), so for a clean run you have to manually
delete the tests' directories.
"""
import argparse
import logging
import os
import shlex
import subprocess
import sys
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import pandas as pd

import apo_holo_structure_stats.settings
from apo_holo_structure_stats.pipeline import (
    chains_for_uniprot_ids,
    download_structures,
    filter_structures,
    make_pairs_lcs,
    run_1struct_analyses,
    run_analyses,
)
from apo_holo_structure_stats.pipeline.filter_structures import StructureProcessor

""" Download unp groups, dl (some) structures, filter structures,

make pairs (but there probably no pairs 
would be made, with first 100 test chains? Then I would have to gather the test data myself, from real make_pairs 
output)
- so I could make the tests with 2 datasets - (subsampled) paper and 300 random chains

run analyses - I would need again some number of pairs -–> I could use those from the paper...
 """

ROOT_OUTPUT_DIR = Path('test_pipeline_output')


logger = logging.getLogger()
logger.level = logging.DEBUG


class TestPipelineBase(TestCase):
    class FileNames:
        CHAINS_WITH_UNP = 'chains.json'
        CHAINS_WITH_UNP_PROCESSED = f'__{CHAINS_WITH_UNP}'
        FILTERED_CHAINS = 'filtered_chains.json'
        PAIRS = 'pairs.json'
        ONE_STRUCT_ANALYSES_DIR = 'one_struct_analyses'  # (this script is after pairs, as the set of pdb structures
        # was greatly reduced and here we run API requests to slow and potentially unstable APIs)
        PAIR_ANALYSES = 'pair_analyses.json'

    """unittests are run in alphabetical order. The script names are in the order, too. I didn't want to make a single
    monolithic test, so that you can rerun tests (of course, if their depedencies=previous in order tests were already
    run)"""

    # like in setup.cfg entry points
    COMMANDS_TO_FUNCTIONS = {
        'ah-chains-uniprot': chains_for_uniprot_ids.main,
        'ah-download-structures': download_structures.main,
        'ah-filter-structures': filter_structures.main,
        'ah-make-pairs': make_pairs_lcs.main,
        'ah-run-1struct-analyses': run_1struct_analyses.main,
        'ah-run-analyses': run_analyses.main,
    }

    def run_command_with_mock_argparse(self, command: str):
        print('running:', command)
        argv = shlex.split(command)

        # patch arguments
        script_fn = self.COMMANDS_TO_FUNCTIONS[argv[0]]
        with patch.object(sys, 'argv', argv):
            script_fn()
            # todo assert no/not too many errors (in stderr)? Test almost always pass, as I catch the errors in the
            #  scripts (for good reasons)

    # def run_command_assert_ok(self, command: str):
    #     print('running:', command)
    #     command = shlex.split(command)
    #
    #     result = subprocess.run(command)
    #     self.assertEqual(0, result.returncode)

    def get_download_structures_input(self) -> pd.DataFrame:
        return pd.read_json(self.FileNames.CHAINS_WITH_UNP)

    def test_download_structures(self):
        # for changing the pdb struct download dir (by default is cwd/pdb_structs) see settings.Settings

        input_df = self.get_download_structures_input()
        input_df.to_json(self.FileNames.CHAINS_WITH_UNP_PROCESSED)

        command = f"""\
ah-download-structures --debug --download_threads 2 "{self.FileNames.CHAINS_WITH_UNP_PROCESSED}"
"""
        self.run_command_with_mock_argparse(command)

    def test_filter_structures(self):
        command = f"""\
ah-filter-structures --debug --disallow_download "{self.FileNames.CHAINS_WITH_UNP_PROCESSED}" "{self.FileNames.FILTERED_CHAINS}"
"""
        self.run_command_with_mock_argparse(command)

    def test_make_pairs(self):
        command = f"""\
ah-make-pairs --debug "{self.FileNames.FILTERED_CHAINS}" "{self.FileNames.PAIRS}"
"""
        self.run_command_with_mock_argparse(command)

    def test_run_1struct_analyses(self):
        command = f"""\
ah-run-1struct-analyses --debug --workers 8 "{self.FileNames.PAIRS}"
"""
        self.run_command_with_mock_argparse(command)

    def test_run_analyses(self):
        command = f"""\
ah-run-analyses --debug --workers 8 "{self.FileNames.PAIRS}"
"""
        self.run_command_with_mock_argparse(command)

    @classmethod
    def setUpClass(cls) -> None:
        # dont delete anything..
        # # prepare clean output dir
        # # (and when tests end, it is not deleted, only when the tests are run again)
        # if self.OUTPUT_DIR.exists():
        #     shutil.rmtree(self.OUTPUT_DIR)

        cls.OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
        cls.original_wd = os.getcwd()
        os.chdir(cls.OUTPUT_DIR)

    @classmethod
    def tearDownClass(cls) -> None:
        os.chdir(cls.original_wd)


class TestPipelinePaperDataset(TestPipelineBase):
    OUTPUT_DIR = ROOT_OUTPUT_DIR / 'paper_dataset'

    def test_chains_for_uniprot_ids(self):
        # the following file does not have uniprotkb ids for the chains,
        # but it will be added by the script
        paper_chains_json_path = Path(self.original_wd) / 'test_data/paper_chains.json'

        command = f"""\
ah-chains-uniprot --debug --chains "{paper_chains_json_path}" "{self.FileNames.CHAINS_WITH_UNP}"
"""
        self.run_command_with_mock_argparse(command)

    # the following definitions are only for the run test button to show up in PyCharm
    def test_download_structures(self):
        super().test_download_structures()

    def test_filter_structures(self):
        super().test_filter_structures()

    def test_make_pairs(self):
        super().test_make_pairs()

    def test_run_1struct_analyses(self):
        super().test_run_1struct_analyses()

    def test_run_analyses(self):
        super().test_run_analyses()


class TestPipelineRandomDataset(TestPipelineBase):
    OUTPUT_DIR = ROOT_OUTPUT_DIR / 'random_dataset'
    SEED = 42
    DATASET_SIZE__CHAINS = 100

    def get_download_structures_input(self) -> pd.DataFrame:
        df = pd.read_json(self.FileNames.CHAINS_WITH_UNP)
        return df.sample(self.DATASET_SIZE__CHAINS, random_state=self.SEED)

    def test_chains_for_uniprot_ids(self):
        command = f"""\
ah-chains-uniprot --debug "{self.FileNames.CHAINS_WITH_UNP}"
"""
        self.run_command_with_mock_argparse(command)

    # the following definitions are only for the run test button to show up in PyCharm
    def test_download_structures(self):
        super().test_download_structures()

    def test_filter_structures(self):
        super().test_filter_structures()

    def test_make_pairs(self):
        super().test_make_pairs()

    # the following tests may fail as in 100 or so random chains there is small chance there will be any pairs
    # and make_pairs won't create an output file if its output were to be empty.
    def test_run_1struct_analyses(self):
        super().test_run_1struct_analyses()

    def test_run_analyses(self):
        super().test_run_analyses()


class TestFilterStructures(TestCase):
    def test_elmi_structure(self):
        """5A1A"""
        with patch('apo_holo_structure_stats.settings.Settings.MIN_STRUCTURE_RESOLUTION', 2.5):
            # has _em_3d_reconstruction.resolution 2.2
            chains_metadata = StructureProcessor().process_structure(Path('test_data/5a1a.cif.gz'))
        self.assertEqual(4, len(chains_metadata))

        with patch('apo_holo_structure_stats.settings.Settings.MIN_STRUCTURE_RESOLUTION', 2.1):
            # has _em_3d_reconstruction.resolution 2.2
            chains_metadata = StructureProcessor().process_structure(Path('test_data/5a1a.cif.gz'))
        self.assertEqual(0, len(chains_metadata))
