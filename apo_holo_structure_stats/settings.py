from pathlib import Path

# STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path('pdb_structs')  # location for Bio.PDB.PDBList that caches the downloaded structures.
STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path('apo_holo_structure_stats/paper_repl/pdb_structs')  # location for Bio.PDB.PDBList that caches the downloaded structures.
# Saves into subdirectories named by the 2 middle letters of a pdb code. Example: 2a3f -> a3/2a3f.cif

MIN_STRUCTURE_RESOLUTION = 2.5  # todo co považuje biopython za resolution

# todo isholo params
