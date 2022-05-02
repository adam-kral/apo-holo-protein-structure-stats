import os
import re
from pathlib import Path

from submit_run_analyses import submit_job, JOB_SCRIPTS_PATH

out_basename = 'filter_output.json'

for job_script_path in JOB_SCRIPTS_PATH.glob('job_script*'):  # zip range for testing
    shard_num = re.findall(r'\d+', job_script_path.name)[0]

    if os.path.exists(Path(out_basename + '.d') / f'{out_basename}{shard_num}'):
        continue
        # successful job (pip got installed and the program ran)

    # print(job_script_path)
    submit_job(job_script_path, 1, '4gb', walltime='8:00:00')


    #
    #
    # # for failed cryptic < 2 files (out gz, json out)
    # # if len(os.listdir(f'output/{shard_num}')) < 2:
    #
    #
    #
    # # less than 3 files (out gz, apo_holo, domains)
    # if len(os.listdir(f'output/{shard_num}')) < 3:
    #
