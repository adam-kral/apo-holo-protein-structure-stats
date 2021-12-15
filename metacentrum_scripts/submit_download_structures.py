import argparse
import os
import subprocess
from pathlib import Path

from run_pipeline import get_shell_template, ShellTemplate
from submit_filter_structures import get_storage_path

base_template = get_shell_template()

init_vars_template = ShellTemplate('''
CODE_DIR=<><code_dir
JOB_OUTPUT_DIR=<><job_output_dir
JOB_OUTPUT_DIR__IN_STORAGE_HOME=<><job_output_dir_in_home
''')

run_script_template = ShellTemplate('''
INPUT_FILE_DIR=../input
mkdir -p "$INPUT_FILE_DIR"
INPUT_FILE=$INPUT_FILE_DIR/<><input_file_name
cp "<><input_file_cp_path" "$INPUT_FILE"

ah-download-structures <><script_opts "$INPUT_FILE"
''')

STORAGE_DIR = Path('.')
JOBS_SCRIPTS_DIR = Path(STORAGE_DIR, 'job_scripts')
JOBS_SCRIPTS_DIR.mkdir(parents=True)

parser = argparse.ArgumentParser()
parser.add_argument('--script_opts', default='--download_threads 100 --debug')
# parser.add_argument('--pdb_dir', default='pdb_structs', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')

parser.add_argument('input_json', type=Path, help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
parser.add_argument('code_dir', type=Path, help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')

args = parser.parse_args()


job_output_dir = get_storage_path(STORAGE_DIR)
init_vars = init_vars_template.substitute(
    code_dir=get_storage_path(args.code_dir),
    job_output_dir=job_output_dir,
    job_output_dir_in_home=job_output_dir.relative_to(os.environ['HOME']),
)

run = run_script_template.substitute(
    input_file_name='chains.json',

    input_file_cp_path=get_storage_path(args.input_json),
    script_opts=args.script_opts,
)

job_script = base_template.substitute(
    init_vars=init_vars,
    run_in_output_dir=run
)

# save script
# could instead pass the script directly to qsub stdin, but this is more explicit
job_script_path = Path(JOBS_SCRIPTS_DIR, f'job_script.sh')
with job_script_path.open('w') as f:
    f.write(job_script)

# submit job
qsub = ['qsub', '-l', 'select=1:ncpus=2:mem=16gb:scratch_local=100gb', '-l', 'walltime=24:00:00', job_script_path]
subprocess.run(qsub)


# download speed can be measured by file size
