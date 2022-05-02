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
INPUT_FILE=$INPUT_FILE_DIR/pairs.json
cp "<><input_file_cp_path" "$INPUT_FILE"

ah-run-1struct-analyses <><script_opts "$INPUT_FILE"
''')

STORAGE_DIR = Path('.')
JOBS_SCRIPTS_DIR = Path(STORAGE_DIR, 'job_scripts')
JOBS_SCRIPTS_DIR.mkdir(parents=True)

parser = argparse.ArgumentParser()
parser.add_argument('--script_opts', default='--debug')

parser.add_argument('input_json', type=Path, help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
parser.add_argument('code_dir', type=Path, help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')

args = parser.parse_args()


job_output_dir = get_storage_path(STORAGE_DIR)

run = run_script_template.substitute(
    input_file_cp_path=get_storage_path(args.input_json),
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
job_script_path = Path(JOBS_SCRIPTS_DIR, f'job_script.sh')
with job_script_path.open('w') as f:
    f.write(job_script)

submit_job(job_script_path, 1, '16gb', scratch_local='10gb')
