import argparse
import json
import logging
import math
import os
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

init_vars_template = ShellTemplate('''
CODE_DIR=<><code_dir
JOB_OUTPUT_DIR=<><job_output_dir
JOB_OUTPUT_DIR__IN_STORAGE_HOME=<><job_output_dir_in_home
SHARD_NUM=<><shard_num
''')

run_script_template = ShellTemplate('''
INPUT_FILE_DIR=../input
mkdir -p "$INPUT_FILE_DIR"

INPUT_FILE=$INPUT_FILE_DIR/<><input_file_name
cp "<><input_shard_dir_cp_path/<><input_file_name" "$INPUT_FILE"

PDB_ARCHIVE=<><{input_shard_pdb_dir}.tar
ssh $(whoami)@storage-brno2.metacentrum.cz "tar -C ~/<><input_shard_dir__in_home -cvf ~/<><input_shard_dir__in_home/$PDB_ARCHIVE <><{input_shard_pdb_dir}"
cp "<><input_shard_dir_cp_path/$PDB_ARCHIVE" "$INPUT_FILE_DIR"
tar -C "$INPUT_FILE_DIR" -xvf "$INPUT_FILE_DIR/$PDB_ARCHIVE"

export STRUCTURE_DOWNLOAD_ROOT_DIRECTORY=$INPUT_FILE_DIR/<><input_shard_pdb_dir

ah-filter-structures <><script_opts "$INPUT_FILE" "<><output_file_name" \
  || { echo >&2 "ah-filter-structures failed (with a code $?)" ;}
''')

STORAGE_DIR = Path('.')
JOBS_SCRIPTS_DIR = Path(STORAGE_DIR, 'job_scripts')


def get_storage_path(relative_path: Path):
    if relative_path.is_absolute():
        return relative_path

    return Path(os.environ['PWD'], relative_path)


def submit_filter_structures(chains, jobs: int, input_shard_base_name: str, output_shard_base_name: str,
                             script_opts: str, code_dir: str, pdb_dir: str):
    shard_size = math.ceil(len(chains) / jobs)
    shards = (chains[i:min(len(chains), i+shard_size)] for i in range(0, len(chains), shard_size))

    # convert 1 -> '01', (pad with zeros to max job_num digit count)
    format_num = lambda n: f'{n:0{len(str(jobs - 1))}}'

    code_dir = Path(code_dir)
    pdb_dir = Path(pdb_dir)
    input_dir = Path(STORAGE_DIR, input_shard_base_name + '.d')
    output_dir = Path(STORAGE_DIR, output_shard_base_name + '.d')

    # todo možná input a output, input->job_scripts, input->json_shards, input->pdb_structs, output -> ...

    # pdb_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    JOBS_SCRIPTS_DIR.mkdir(parents=True)

    for job_num, shard in enumerate(shards):
        # save shard
        input_shard_path = Path(input_dir, input_shard_base_name + format_num(job_num))
        with input_shard_path.open('w') as f:
            json.dump(shard, f)

        # prepare input structure directory for a job (input shard)
        # the job will tar it (by ssh'ing to the storage server) and copy it to its SCRATCHDIR
        # OR do this in the job? os.link should be fast and frontend has the storage as default home dir
        # what is resource intensive is the qsub part? (Or that just waits and scheduling is done elswhere?)
        input_shard_pdb_subdir = Path(input_dir, 'pdb_structs' + format_num(job_num))
        input_shard_pdb_subdir.mkdir()
        for pdb_code in set(chain['pdb_code'] for chain in shard):
            fname = f'{pdb_code}.cif.gz'
            os.link(pdb_dir / fname, input_shard_pdb_subdir / fname)

        # render job script (for the execution host)
        # render the `run_in_output_dir` block
        output_shard_path = Path(output_dir, output_shard_base_name + format_num(job_num))

        run = run_script_template.substitute(
            output_file_name=output_shard_path.name,

            input_file_name=input_shard_path.name,
            input_shard_dir_cp_path=get_storage_path(input_dir),
            input_shard_dir__in_home=get_storage_path(input_dir).relative_to(os.environ['HOME']),

            input_shard_pdb_dir=input_shard_pdb_subdir.name,

            script_opts=script_opts,
        )
        job_output_dir = get_storage_path(output_shard_path.parent)
        init_vars = init_vars_template.substitute(
            code_dir=get_storage_path(code_dir),
            job_output_dir=job_output_dir,
            job_output_dir_in_home=job_output_dir.relative_to(os.environ['HOME']),
            shard_num=job_num,
        )
        job_script = base_template.substitute(
            init_vars=init_vars,
            run_in_output_dir=run
        )

        # save script
        # could instead pass the script directly to qsub stdin, but this is more explicit
        job_script_path = Path(JOBS_SCRIPTS_DIR, f'job_script{format_num(job_num)}.sh')
        with job_script_path.open('w') as f:
            f.write(job_script)

        # submit job
        # passed_env_vars = [
        #     f'CODE_DIR={get_storage_path(code_dir)}',
        #     f'JOB_OUTPUT_DIR={get_storage_path(output_shard_path.parent)}',
        #     f'STRUCTURE_DOWNLOAD_ROOT_DIRECTORY={get_storage_path(pdb_dir)}',
        # ]
        # passed_env_vars = ','.join(passed_env_vars)  # no <space> between ',' Shouldn't have used it in the manpage, if it doesn't work!!
        # large walltime - sometimes copying takes long?? Wtf, why?
        # qsub = ['qsub', '-l', 'select=1:ncpus=4:mem=8gb:scratch_local=10gb', '-l', 'walltime=4:00:00', job_script_path]
        # logger.info('running qsub: ' + ' '.join(map(str, qsub)))
        # break  # for testing, run only one job


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--jobs', type=int, default=1, help='number of jobs')
    # todo add reasonable default
    parser.add_argument('--script_opts', default='--disallow_download --workers 4 --debug')
    # parser.add_argument('--pdb_dir', default='pdb_structs', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')

    parser.add_argument('input_json', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('pdb_dir', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('code_dir', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    # todo following two args might be optional
    parser.add_argument('input_shard_base_name', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('output_shard_base_name', help='output filename for the json list of pdb_codes that passed the filter. Paths to mmcif files are relative to the working directory.')
    add_loglevel_args(parser)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)  # bohužel musim specifikovat i tohle, protoze takhle to s __name__ funguje...
    logging.basicConfig()

    with open(args.input_json) as f:
        chains = json.load(f)

    submit_filter_structures(chains, args.jobs, args.input_shard_base_name, args.output_shard_base_name,
                             args.script_opts, args.code_dir, args.pdb_dir)
