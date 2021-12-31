import subprocess
from pathlib import Path

JOB_SCRIPTS_PATH = Path('job_scripts')

for job_script_path in JOB_SCRIPTS_PATH.glob('job_script*'):
        qsub = ['qsub', '-l', 'select=1:ncpus=4:mem=12gb:scratch_local=10gb', '-l', 'walltime=12:00:00', job_script_path]
        print('running qsub: ' + ' '.join(map(str, qsub)))
        subprocess.run(qsub)
        # break  # for testing, run only one job
