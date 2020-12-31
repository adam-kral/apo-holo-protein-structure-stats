# apo-holo-protein-structure-stats

## Install
In repository root:
- Create a new virtual environment (recommended) named 'venv' with Python 3: `python3 -m venv venv` (Python 3 has a built-in `venv` module for that, alternatively you can use virtualenv: `virtualenv -p python3 venv`)
- Activate the venv: `source venv/bin/activate` on Unix or `venv\Scripts\activate.bat` on Windows
- Install requirements `pip install -r requirements.txt`
- todo for easier script running, extract runnable scripts to the repository root/maybe setup.py

## Run (and pipeline description)
Run following commands in your virtual environment and in repository root.
(Script are run as Python modules within the package (`-m` flag). That's because they import from the package. In future I could extract
them to the repository root (just the  `if __name__=='__main__'` code), so they could be run just by their name/path.)


### Filter structures
Filters structures suitable for analysis. 

Downloads the structures with codes given in `pdb_codes_or_directory` argument, or traverses
 the directory tree at `pdb_codes_or_directory` if `-d` option is present (all files are assumed to be mmcif files). It loads the structures
 and decides if structure meets criteria for resolution, single-chainedness, etc. In the json output file, there are only structures which
 passed the filter.

todo - usage is rather `python -m apo_holo_structure_stats.pipeline.filter_structures`
```
usage: filter_structures.py [-h] [-d] pdb_codes_or_directory output_file

positional arguments:
  pdb_codes_or_directory
                        comma-delimited list of pdb_codes, or if `-d` option
                        is present, a directory with mmcif files.
  output_file           output filename for the json list of pdb_codes that
                        passed the filter

optional arguments:
  -h, --help            show this help message and exit
  -d, --is-directory    now pdb_codes_or_directory is a path to a directory
                        with mmcif files. Whole tree structure is inspected,
                        all files are assumed to be mmcifs.
```

Example command:

`python -m apo_holo_structure_stats.pipeline.filter_structures 3mqd,3lrf,4jv3,3u0e,3u0f o_filter.json`  
or  
`python -m apo_holo_structure_stats.pipeline.filter_structures -d my_directory_with_mmcifs o_filter.json`

### Annotate with ligand presence bool (is_holo)

```
usage: is_holo.py [-h] structures_json output_file

positional arguments:
  structures_json  annotate the list of structures with is_holo bool. File
                   needs to contain list of objects with pdb_code and
                   main_chain_id and path to the structure
  output_file      writes input json annotated with boolean "is holo"

optional arguments:
  -h, --help       show this help message and exit
```

### Annotate with isoform information

```
usage: isoform.py [-h] structures_json output_file

positional arguments:
  structures_json  annotate the list of structures with isoform data. File
                   needs to contain list of objects with pdb_code and
                   main_chain_id
  output_file      writes input json annotated with isoform uniprot id

optional arguments:
  -h, --help       show this help message and exit

```

### Run analyses
```
usage: run_analyses.py [-h] [--isoform ISOFORM] structures_json output_file

positional arguments:
  structures_json    list of structures {pdb_code: , path: , isoform_id: ,
                     is_holo: bool, ?main_chain_id: }
  output_file        dumped results of analyses

optional arguments:
  -h, --help         show this help message and exit
  --isoform ISOFORM  process only structures with main chain of that isoform
```

### Try out the whole pipeline
```shell script
python -m apo_holo_structure_stats.pipeline.filter_structures 3mqd,3lrf,4jv3,3u0e,3u0f filter.json
python -m apo_holo_structure_stats.pipeline.is_holo filter.json is_holo.json
python -m apo_holo_structure_stats.pipeline.isoform is_holo.json structures.json
python -m apo_holo_structure_stats.pipeline.run_analyses structures.json analyses_output.json
```
