import subprocess
from pathlib import Path


def submit_job(script_path):
    qsub = ['qsub', '-l', 'select=1:ncpus=4:mem=12gb:scratch_local=10gb', '-l', 'walltime=12:00:00', script_path]
    print('running qsub: ' + ' '.join(map(str, qsub)))
    subprocess.run(qsub)


JOB_SCRIPTS_PATH = Path('job_scripts')


if __name__ == '__main__':
    for job_script_path in JOB_SCRIPTS_PATH.glob('job_script*'):
        submit_job(job_script_path)
        # break  # for testing, run only one job
