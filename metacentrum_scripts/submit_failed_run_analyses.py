import os
import re

from submit_run_analyses import submit_job, JOB_SCRIPTS_PATH


for job_script_path in JOB_SCRIPTS_PATH.glob('job_script*'):
    shard_num = re.findall(r'\d+', job_script_path.name)[0]
    # less than 3 files (out gz, apo_holo, domains)
    if len(os.listdir(f'output/{shard_num}')) < 3:
        submit_job(job_script_path)
