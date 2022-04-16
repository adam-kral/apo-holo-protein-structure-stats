import argparse
import itertools
import json
import logging
import math
import os
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


# todo not [ ] output or [x] input file, but always copy dir
# - cp -r is ok.. (
base_template = get_shell_template()

init_vars_template = ShellTemplate('''
CODE_DIR=<><code_dir
JOB_OUTPUT_DIR=<><job_output_dir
JOB_OUTPUT_DIR__IN_STORAGE_HOME=<><job_output_dir_in_home
SHARD_NUM=<><shard_num
''')


#  todo in storage home important!
run_script_template = ShellTemplate('''
INPUT_FILE_DIR=../input
mkdir -p "$INPUT_FILE_DIR"

cp -r <><input_shard_dir_cp_path/* "$INPUT_FILE_DIR"

INPUT_ARCHIVE=input_archive.tar
ssh storage-brno2.metacentrum.cz "tar -C ~/<><input_shard_dir__in_home -cvf ~/<><input_shard_dir__in_home/$INPUT_ARCHIVE ."
cp "<><input_shard_dir_cp_path/$INPUT_ARCHIVE" "$INPUT_FILE_DIR"
tar -C "$INPUT_FILE_DIR" -xvf "$INPUT_FILE_DIR/$INPUT_ARCHIVE"

export STRUCTURE_DOWNLOAD_ROOT_DIRECTORY=$INPUT_FILE_DIR/pdb_structs

# e.g. ah-run-analyses for script_name
<><script_name <><script_opts --opt_input_dir "$INPUT_FILE_DIR" "$INPUT_FILE_DIR/<><input_file_name" \
  || { echo >&2 "<><script_name failed (with a code $?)" ;}
''')

STORAGE_DIR = Path('.')
JOBS_SCRIPTS_DIR = Path(STORAGE_DIR, 'job_scripts')


def get_storage_path(relative_path: Path):
    if relative_path.is_absolute():
        return relative_path

    try:
        # if this script is run in a job
        storage_base = os.environ['PBS_O_WORKDIR']
    except KeyError:
        # this script is run on frontend (not when this scripts runs in a SCRATCHDIR on a computational node)
        storage_base = os.environ['PWD']

    return Path(storage_base, relative_path)


def get_submission_home() -> Path:
    try:
        return Path(os.environ['PBS_O_HOME'])
    except KeyError:
        # this script is run on frontend (not when this scripts runs in a SCRATCHDIR on a computational node)
        return Path(os.environ['HOME'])



def submit_run_analyses(pairs, jobs: int, input_shard_base_name: str, script_opts: str, code_dir: str, pdb_dir: str,
                        qsub_oe_path: Path, other_inputs_dir: str, script_name: str):

    shard_size = math.ceil(len(pairs) / jobs)
    shards = (pairs[i:min(len(pairs), i+shard_size)] for i in range(0, len(pairs), shard_size))

    # convert 1 -> '01', (pad with zeros to max job_num digit count)
    format_num = lambda n: f'{n:0{len(str(jobs - 1))}}'

    code_dir = Path(code_dir)
    pdb_dir = Path(pdb_dir)
    input_dir = Path(STORAGE_DIR, 'input')
    output_dir = Path(STORAGE_DIR, 'output')

    # EDIT NE možná input a output, input->job_scripts, input->json_shards, input->pdb_structs, output -> ...
    # todo input relative to storage
    # pdb_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True)
    output_dir.mkdir(parents=True)
    JOBS_SCRIPTS_DIR.mkdir(parents=True)

    for job_num, shard in enumerate(shards):
        # prepare input structure directory for a job (input shard)
        # the job will tar it (by ssh'ing to the storage server) and copy it to its SCRATCHDIR
        # OR do this in the job? os.link should be fast and frontend has the storage as default home dir
        # what is resource intensive is the qsub part? (Or that just waits and scheduling is done elswhere?)
        output_shard_dir = Path(output_dir, format_num(job_num))
        output_shard_dir.mkdir()

        input_shard_dir = Path(input_dir, format_num(job_num))
        input_shard_dir.mkdir()

        # save shard
        input_shard_path = Path(input_shard_dir, input_shard_base_name)
        with input_shard_path.open('w') as f:
            json.dump(shard, f)

        # for os link to work, we need to link within the same filesystem, so do it
        input_shard_pdb_subdir = get_storage_path(Path(input_shard_dir, 'pdb_structs'))
        input_shard_pdb_subdir.mkdir(parents=True)

        pdb_codes = set(itertools.chain(*([pair['pdb_code_apo'], pair['pdb_code_holo']] for pair in shard)))
        for pdb_code in pdb_codes:
            fname = f'{pdb_code}.cif.gz'
            os.link(pdb_dir / fname, input_shard_pdb_subdir / fname)

        for file in Path(other_inputs_dir).iterdir():
            if not file.is_file():
                continue

            os.link(file, get_storage_path(input_shard_dir) / file.name)

        # render job script (for the execution host)
        # render the `run_in_output_dir` block

        run = run_script_template.substitute(
            input_file_name=input_shard_path.name,
            input_shard_dir_cp_path=get_storage_path(input_shard_dir),
            input_shard_dir__in_home=get_storage_path(input_shard_dir).relative_to(get_submission_home()),

            input_shard_pdb_dir=input_shard_pdb_subdir.name,

            script_opts=script_opts,
            script_name=script_name,
        )
        job_output_dir = get_storage_path(output_shard_dir)
        init_vars = init_vars_template.substitute(
            code_dir=get_storage_path(code_dir),
            job_output_dir=job_output_dir,
            job_output_dir_in_home=job_output_dir.relative_to(get_submission_home()),
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

        logger.info(f'job prepared for shard {job_num} / {jobs}')

        # submit job
        # passed_env_vars = [
        #     f'CODE_DIR={get_storage_path(code_dir)}',
        #     f'JOB_OUTPUT_DIR={get_storage_path(output_shard_path.parent)}',
        #     f'STRUCTURE_DOWNLOAD_ROOT_DIRECTORY={get_storage_path(pdb_dir)}',
        # ]
        # passed_env_vars = ','.join(passed_env_vars)  # no <space> between ',' Shouldn't have used it in the manpage, if it doesn't work!!
        # large walltime - sometimes copying takes long?? Wtf, why?



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--script_name', default='ah-run-analyses')
    parser.add_argument('--jobs', type=int, default=1, help='number of jobs')
    # todo add reasonable default
    parser.add_argument('--script_opts', default='--workers 4 --debug')
    parser.add_argument('--qsub_oe_path', type=Path, default=Path())
    # parser.add_argument('--pdb_dir', default='pdb_structs', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')

    parser.add_argument('input_dir', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('input_json', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('pdb_dir', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    parser.add_argument('code_dir', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    # todo following two args might be optional
    parser.add_argument('input_shard_base_name', help='comma-delimited list of pdb_codes, or if `-d` option is present, a directory with mmcif files.')
    add_loglevel_args(parser)

    args = parser.parse_args()
    logger.setLevel(args.loglevel)  # bohužel musim specifikovat i tohle, protoze takhle to s __name__ funguje...
    logging.basicConfig()

    with open(args.input_json) as f:
        pairs = json.load(f)

    # [done] ttodo zapomnel jsem vyfiltrovat pary s mismatches, takze mam joby, co bezej strasne dlouho nebo joby, co jsou rychly...
    #  ty se filtrujou az v run analysis, nektery joby mely 5700/7500 a jiny treba 2?
    pairs = [pair for pair in pairs if pair['lcs_result']['mismatches'] == 0]

    submit_run_analyses(pairs, args.jobs, args.input_shard_base_name, args.script_opts, args.code_dir, args.pdb_dir,
                        args.qsub_oe_path, args.input_dir, args.script_name)
