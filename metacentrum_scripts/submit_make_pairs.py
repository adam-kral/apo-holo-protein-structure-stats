import argparse
import os
from pathlib import Path

from submit_run_analyses import submit_job
from run_pipeline import get_base_run_script_template, ShellTemplate
from submit_filter_structures import get_storage_path

base_template = get_base_run_script_template()

run_script_template = ShellTemplate('''
INPUT_FILE_DIR=../input
mkdir -p "$INPUT_FILE_DIR"

# copy there the dir of filter structures output
ARCHIVE=input_make_pairs.tar.gz
# gzip is not there, so cvf, not czvf
ssh storage-brno2.metacentrum.cz "tar -C ~/<><input_dir_path__in_storage_home -cvf ~/<><input_dir_path__in_storage_home/$ARCHIVE ."
cp "<><input_dir_path/$ARCHIVE" "$INPUT_FILE_DIR"
mkdir -p "$INPUT_FILE_DIR/filtered"
tar -C "$INPUT_FILE_DIR/filtered" -xvf "$INPUT_FILE_DIR/$ARCHIVE"

ah-make-pairs <><script_opts "$INPUT_FILE_DIR"'/filtered/*.json*' pairs.json
''')

STORAGE_DIR = Path('.')
JOBS_SCRIPTS_DIR = Path(STORAGE_DIR, 'job_scripts')
JOBS_SCRIPTS_DIR.mkdir(parents=True)

parser = argparse.ArgumentParser()
parser.add_argument('--script_opts', default='--debug --workers 8')
# parser.add_argument('--pdb_dir', default='pdb_structs', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')

parser.add_argument('filter_structures_output_dir', type=Path, help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
parser.add_argument('chains_json', type=Path, help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
parser.add_argument('code_dir', type=Path, help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')

args = parser.parse_args()


job_output_dir = get_storage_path(STORAGE_DIR)

run = run_script_template.substitute(
    input_file_cp_path=get_storage_path(args.chains_json),
    input_dir_path=get_storage_path(args.filter_structures_output_dir),
    input_dir_path__in_storage_home=get_storage_path(args.filter_structures_output_dir).relative_to(os.environ['HOME']),

    script_opts=args.script_opts,

)

job_script = base_template.substitute(
    code_dir=get_storage_path(args.code_dir),
    job_output_dir=job_output_dir,
    job_output_dir_in_home=job_output_dir.relative_to(os.environ['HOME']),
    shard_num='',  # N/A

    run_in_output_dir=run,
)

# save script
# could instead pass the script directly to qsub stdin, but this allows one to inspect it
job_script_path = Path(JOBS_SCRIPTS_DIR, f'job_script.sh')
with job_script_path.open('w') as f:
    f.write(job_script)

submit_job(job_script_path, 8, '24gb', scratch_local='10gb')  # took an hour without multiprocessing or pypy on 1 cpu

# download speed can be measured by file size
