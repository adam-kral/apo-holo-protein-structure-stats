from pathlib import Path

# location for Bio.PDB.PDBList that caches the downloaded structures. Saves into subdirectories named by the 2 middle letters of a pdb code
#   2a3f -> a3/2a3f.cif
STRUCTURE_DOWNLOAD_ROOT_DIRECTORY = Path.cwd() / 'pdb_structs'

# todo co považuje biopython za resolution
MIN_STRUCTURE_RESOLUTION = 2.5
