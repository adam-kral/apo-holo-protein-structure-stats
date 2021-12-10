# move to execution-local directory
# first check if set (is done automatically)
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }
cd "$SCRATCHDIR"

# these two variables I pass via qsub -v option
test -n "$CODE_DIR" || { echo >&2 "Variable CODE_DIR is not set!"; exit 1; }
test -n "$JOB_OUTPUT_DIR" || { echo >&2 "Variable JOB_OUTPUT_DIR is not set!"; exit 1; }

# install requirements
module add python/3.8.0-gcc-rab6t
python -m venv venv
source venv/bin/activate
#python -m ensurepip --upgrade  # upgrade uplne nefunguval?  # edit: vypada to, ze ve venvu se pip vytvori
pip install --upgrade pip  # proto

# copy project sources and install the project in "src" dir
cp -r "$CODE_DIR" .
cd "$(basename "$CODE_DIR")"
pip install -r requirements.txt
# this command will install the apo/holo project in site-packages,
# so we can easily run project's scripts directly
pip install -e .
cd "$SCRATCHDIR"

# move to directory for script output
output_dir=output
mkdir "$output_dir"
cd "$output_dir"

# TEMPLATE
# For example:
# - (optionally) copy input files for the program here
# - run a script with some arguments
><>run_in_output_dir

# copy output to (non-execution host) permanent storage
cp -r . "$JOB_OUTPUT_DIR" || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# clean the SCRATCH directory
clean_scratch
