from pathlib import Path
import os

# location for Bio.PDB.PDBList that caches the downloaded structures.
# Saves into subdirectories named by the 2 middle letters of a pdb code. Example: 2a3f -> a3/2a3f.cif
STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path(os.environ.get('STRUCTURE_DOWNLOAD_ROOT_DIRECTORY', 'pdb_structs'))


# STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path(Path(__file__).parent, 'paper_repl', 'pdb_structs')

MIN_STRUCTURE_RESOLUTION = 2.5  # todo co považuje biopython za resolution

# todo isholo params
