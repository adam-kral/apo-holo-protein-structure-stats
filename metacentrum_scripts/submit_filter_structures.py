import argparse
import json
import logging
import math
import subprocess
from pathlib import Path

from run_pipeline import get_shell_template, ShellTemplate

logger = logging.getLogger(__name__)

def add_loglevel_args(parser):
    parser.add_argument(
        '--debug',
        action="store_const", dest="loglevel", const=logging.DEBUG,
        default=logging.WARNING,
    )
    parser.add_argument(
        '-v', '--verbose',
        action="store_const", dest="loglevel", const=logging.INFO,
    )



base_template = get_shell_template()
run_script_template = ShellTemplate('''
INPUT_FILE_DIR=../input
mkdir -p "$INPUT_FILE_DIR"
INPUT_FILE=$INPUT_FILE_DIR/><>input_file_name
cp "><>input_file_cp_path" "$INPUT_FILE"

ah-filter-structures ><>script_opts "$INPUT_FILE" "><>output_file_name"
''')

STORAGE_DIR = Path('.')
JOBS_SCRIPTS_DIR = Path(STORAGE_DIR, 'job_scripts')


def submit_filter_structures(chains, jobs: int, input_shard_base_name: str, output_shard_base_name: str,
                             script_opts: str, code_dir: str):
    shard_size = math.ceil(len(chains) / jobs)
    shards = (chains[i:min(len(chains), i+shard_size)] for i in range(0, len(chains), shard_size))

    # convert 1 -> '01', (pad with zeros to max job_num digit count)
    format_num = lambda n: f'{n:0{len(str(jobs - 1))}}'

    code_dir = Path(code_dir)
    input_dir = Path(STORAGE_DIR, input_shard_base_name + '.d')
    output_dir = Path(STORAGE_DIR, output_shard_base_name + '.d')

    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    JOBS_SCRIPTS_DIR.mkdir(parents=True)

    for job_num, shard in enumerate(shards):
        # save shard
        input_shard_path = Path(input_dir, input_shard_base_name + format_num(job_num))
        with input_shard_path.open('w') as f:
            json.dump(shard, f)

        # render job script (for the execution host)
        # render the `run_in_output_dir` block
        output_shard_path = Path(output_dir, output_shard_base_name + format_num(job_num))

        run = run_script_template.substitute(
            output_file_name=output_shard_path.name,
            input_file_name=input_shard_path.name,

            input_file_cp_path=input_shard_path.absolute(),
            script_opts=script_opts,
        )
        job_script = base_template.substitute(run_in_output_dir=run)

        # save script
        # could instead pass the script directly to qsub stdin, but this is more explicit
        job_script_path = Path(JOBS_SCRIPTS_DIR, f'job_script{format_num(job_num)}.sh')
        with job_script_path.open('w') as f:
            f.write(job_script)

        # submit job
        passed_env_vars = [
            f'CODE_DIR={code_dir.absolute()}',
            f'JOB_OUTPUT_DIR={output_shard_path.parent.absolute()}',
        ]
        passed_env_vars = ','.join(passed_env_vars)  # no <space> between ',' Shouldn't have used it in the manpage, if it doesn't work!!
        qsub = ['qsub', '-l', 'select=1:ncpus=4:mem=2gb:scratch_local=10gb', '-l', 'walltime=2:00:00', '-v', passed_env_vars, job_script_path]
        logger.info('running qsub: ' + ' '.join(map(str, qsub)))
        subprocess.run(qsub)
        break  # for testing, run only one job


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--jobs', type=int, default=1, help='number of jobs')
    # todo add reasonable default
    parser.add_argument('--script_opts', default='--download_threads 10 --workers 4 --debug')

    parser.add_argument('code_dir', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    # todo following two args might be optional
    parser.add_argument('input_shard_base_name', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('output_shard_base_name', help='output filename for the json list of pdb_codes that passed the filter. Paths to mmcif files are relative to the working directory.')
    parser.add_argument('input_json', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    add_loglevel_args(parser)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)  # bohužel musim specifikovat i tohle, protoze takhle to s __name__ funguje...
    logging.basicConfig()

    with open(args.input_json) as f:
        chains = json.load(f)

    submit_filter_structures(chains, args.jobs, args.input_shard_base_name, args.output_shard_base_name,
                             args.script_opts, args.code_dir)
