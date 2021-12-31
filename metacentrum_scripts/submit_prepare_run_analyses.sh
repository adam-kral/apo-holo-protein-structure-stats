STORAGE_HOME=/storage/brno2/home/$(whoami)
AH_OUTPUT_DIR=$STORAGE_HOME/ah_output/1212

METACENTRUM_SCRIPT_DIR=$STORAGE_HOME/apo-holo-protein-structure-stats/metacentrum_scripts
JOB_OUTPUT_PATH=$PBS_O_WORKDIR

cd "$SCRATCHDIR"

module add python/3.8.0-gcc-rab6t

JOBS=400
python "$METACENTRUM_SCRIPT_DIR/prepare_run_analyses.py" --debug --qsub_oe_path "$JOB_OUTPUT_PATH" --jobs $JOBS \
  $STORAGE_HOME/apo-holo-protein-structure-stats/run_analyses_input_dir \
  $STORAGE_HOME/apo-holo-protein-structure-stats/o_make_pairs_lcs $AH_OUTPUT_DIR/download/pdb_structs/ \
  $STORAGE_HOME/apo-holo-protein-structure-stats/ pairs.json_shard


# copy output job scripts and json shards
# no, set qsub output location !!
cp -r . "$JOB_OUTPUT_PATH"  || { echo >&2 "Result file(s) copying failed (with a code $?)"; exit 4; }
clean_scratch
