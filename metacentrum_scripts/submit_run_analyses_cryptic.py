import subprocess
from pathlib import Path


def submit_job(script_path, ncpus, mem, walltime='12:00:00', scratch_local='10gb'):
    qsub = ['qsub', '-l', f'select=1:ncpus={ncpus}:mem={mem}:scratch_local={scratch_local}', '-l', f'walltime={walltime}', script_path]
    print('running qsub: ' + ' '.join(map(str, qsub)))
    subprocess.run(qsub)


JOB_SCRIPTS_PATH = Path('job_scripts')


if __name__ == '__main__':
    # for job_script_path, _ in zip(JOB_SCRIPTS_PATH.glob('job_script*'), range(10)):  # zip range for testing
    for job_script_path in JOB_SCRIPTS_PATH.glob('job_script*'):  # zip range for testing
        submit_job(job_script_path, 1, '8gb')
        # break  # for testing, run only one job
