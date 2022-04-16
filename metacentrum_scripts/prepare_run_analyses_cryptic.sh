STORAGE_HOME=/storage/brno2/home/$(whoami)
#SCRIPT_NAME=ah_run_analyses
SCRIPT_NAME=../apo-holo-protein-structure-stats/cryptic_binding_sites/run_analyses.py
AH_OUTPUT_DIR=$STORAGE_HOME/cryptic_output  #zmenit + code_dir nechat, jenom tam prekopirovat cryptic binding sites
PDB_DIR=$STORAGE_HOME/ah_output/1212/download/pdb_structs/
# + zmenit nazev skriptu (celou path nebo pridat do setup.cfg)

METACENTRUM_SCRIPT_DIR=$STORAGE_HOME/apo-holo-protein-structure-stats/metacentrum_scripts
JOB_OUTPUT_PATH=$PBS_O_WORKDIR

cd "$SCRATCHDIR"

module add python/3.8.0-gcc-rab6t

JOBS=3000
python "$METACENTRUM_SCRIPT_DIR/prepare_run_analyses.py" --debug --qsub_oe_path "$JOB_OUTPUT_PATH" --jobs $JOBS \
  --script_name "$SCRIPT_NAME" \
  $STORAGE_HOME/apo-holo-protein-structure-stats/run_analyses_input_dir \
  $STORAGE_HOME/apo-holo-protein-structure-stats/o_make_pairs_lcs $PDB_DIR \
  $STORAGE_HOME/apo-holo-protein-structure-stats/ pairs.json_shard


# copy output job scripts and json shards
# no, set qsub output location !!
cp -r . "$JOB_OUTPUT_PATH"  || { echo >&2 "Result file(s) copying failed (with a code $?)"; exit 4; }
clean_scratch
