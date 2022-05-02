#todo why is this named 'submit' when it does not submit qsub tasks
STORAGE_HOME=/storage/brno2/home/$(whoami)
SCRIPT_NAME=ah-run-analyses

METACENTRUM_SCRIPT_DIR=$STORAGE_HOME/apo-holo-protein-structure-stats/metacentrum_scripts
JOB_OUTPUT_PATH=$PBS_O_WORKDIR

cd "$SCRATCHDIR"

module add python/3.8.0-gcc-rab6t

JOBS=1500  # limit is 10 000 scheduled at once
python "$METACENTRUM_SCRIPT_DIR/prepare_run_analyses.py" --debug --qsub_oe_path "$JOB_OUTPUT_PATH" --jobs $JOBS \
  --script_name "$SCRIPT_NAME" \
  "$JOB_OUTPUT_PATH/../one_struct_analyses" \
  "$JOB_OUTPUT_PATH/../make_pairs/pairs.json" "$JOB_OUTPUT_PATH/../download/pdb_structs/" \
  $STORAGE_HOME/apo-holo-protein-structure-stats/ pairs.json_shard


# copy output job scripts and json shards
# no, set qsub output location !!
cp -r . "$JOB_OUTPUT_PATH"  || { echo >&2 "Result file(s) copying failed (with a code $?)"; exit 4; }
clean_scratch
