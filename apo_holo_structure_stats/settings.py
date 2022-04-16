import dataclasses
from pathlib import Path
import os

# location for Bio.PDB.PDBList that caches the downloaded structures.
# Saves into subdirectories named by the 2 middle letters of a pdb code. Example: 2a3f -> a3/2a3f.cif
STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path(os.environ.get('STRUCTURE_DOWNLOAD_ROOT_DIRECTORY', 'pdb_structs'))

# todo for testing, probably create an env var.. IS_PRODUCTION=1 or something like that
# STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path(Path(__file__).parent, 'paper_repl', 'pdb_structs')

MIN_STRUCTURE_RESOLUTION = 2.5  # todo co považuje biopython za resolution


class LigandSpec:
    MIN_NON_H_ATOMS = 6
    MIN_RESIDUES_WITHIN_LIGAND = 6
    MIN_RESIDUES_WITHIN_LIGAND__RADIUS = 4.5
