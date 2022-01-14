# move to execution-local directory
# first check if set (is done automatically by pbs)
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
cd "$SCRATCHDIR"

<><init_vars

# these two variables I pass via qsub -v option
test -n "$CODE_DIR" || { echo >&2 "Variable CODE_DIR is not set!"; exit 1; }
#test -n "$SHARD_NUM" || { echo >&2 "Variable SHARD_NUM is not set!"; exit 1; }
test -n "$JOB_OUTPUT_DIR" || { echo >&2 "Variable JOB_OUTPUT_DIR is not set!"; exit 1; }
test -n "$JOB_OUTPUT_DIR__IN_STORAGE_HOME" || { echo >&2 "Variable JOB_OUTPUT_DIR__IN_STORAGE_HOME is not set!"; exit 1; }

# install requirements
# python 3.8.0nefungoval
module add python/3.8.0-gcc-rab6t
python -m venv venv
source venv/bin/activate

#module add anaconda3-2019.10
#export CONDA_ALWAYS_YES="true"
#conda update conda  # pujde?
#conda init bash
#exec "$SHELL"  # nefunguje, spadne job
#source ~/.bashrc  # pouzij teda tohle
#conda create -p venv python=3.8
#source activate "$SCRATCHDIR/venv"
python --version

pip install --upgrade pip  # proto

# copy project sources and install the project in "src" dir
cp -r "$CODE_DIR" .
cd "$(basename "$CODE_DIR")"
pip install -r requirements.txt
# this command will install the apo/holo project in site-packages,
# so we can easily run project's scripts directly
pip install -e .
cd "$SCRATCHDIR"

# create TMPDIR in SCRATCH, otherwise outside-scratch quota is (unfortunately sometimes) exceeded (during pip install)
TMP=$SCRATCHDIR/tmp
mkdir "$TMP"
export TMPDIR=$TMP

# move to directory for script output
output_dir=output
mkdir "$output_dir"
cd "$output_dir"

# TEMPLATE
# For example:
# - (optionally) copy input files for the program here
# - run a script with some arguments
<><run_in_output_dir

# copy output to (non-execution host) permanent storage
ARCHIVE=output${SHARD_NUM}.tar.gz
tar -czvf "$ARCHIVE" .
cp "$ARCHIVE" "$JOB_OUTPUT_DIR"  || { echo >&2 "Result file(s) copying failed (with a code $?)"; exit 4; }

# see hardcoded locations at
# https://wiki.metacentrum.cz/wiki/Working_with_data#Data_transfer_between_storages_and_PC.2C_principles
ssh storage-brno2.metacentrum.cz "tar -C ~/$JOB_OUTPUT_DIR__IN_STORAGE_HOME -xvf ~/$JOB_OUTPUT_DIR__IN_STORAGE_HOME/$ARCHIVE" \
  || { echo >&2 "Unpacking archive at dest failed (with a code $?)"; exit 4; }

#cp -r . "$JOB_OUTPUT_DIR" || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# clean the SCRATCH directory
clean_scratch
