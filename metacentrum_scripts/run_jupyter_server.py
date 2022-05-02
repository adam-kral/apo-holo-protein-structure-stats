#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

from prepare_run_analyses import get_storage_path, get_submission_home
from run_pipeline import get_base_run_script_template, ShellTemplate

base_template = get_base_run_script_template()

init_vars_template = ShellTemplate('''
CODE_DIR=<><code_dir
JOB_OUTPUT_DIR=<><job_output_dir
JOB_OUTPUT_DIR__IN_STORAGE_HOME=<><job_output_dir_in_home
''')


run_script_template = ShellTemplate('''
pip install notebook 
# unfortunately remote ntbs don't show ntb variables in pycharm https://youtrack.jetbrains.com/issue/PY-43337

jupyter notebook --port 42007 ../apo-holo-protein-structure-stats/apo_holo_structure_stats/results/results.ipynb \
  || { echo >&2 "jupyter failed (with a code $?)" ;}
''')

# todo ale nebudu znat token!, ssh to comp node and inspect ouput, faster than somehow setting passwd

parser = argparse.ArgumentParser()
parser.add_argument('code_dir', type=Path, help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
args = parser.parse_args()

# cwd in submission host is output dir
job_output_dir = get_storage_path(Path())

init_vars = init_vars_template.substitute(
    code_dir=get_storage_path(args.code_dir),
    job_output_dir=job_output_dir,
    job_output_dir_in_home=job_output_dir.relative_to(get_submission_home()),
)
job_script = base_template.substitute(
    init_vars=init_vars,
    run_in_output_dir=run_script_template,
)
job_script_path = Path('job_script_jupyter.sh   ')

with job_script_path.open('w') as f:
    f.write(job_script)

qsub = ['qsub', '-l', 'select=1:ncpus=3:mem=32gb:scratch_local=10gb', '-l', 'walltime=24:00:00', job_script_path]
print('running qsub: ' + ' '.join(map(str, qsub)))
subprocess.run(qsub)
