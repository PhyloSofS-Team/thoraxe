#!/bin/bash
# Worker of thoraxe/thoraxe_array.py: runs ThorAxe on a share of the genes of OUT_DIR/run/genes.txt.
#
# Do not submit it by hand: thoraxe_array.py prepares OUT_DIR/run/ (code snapshot, aligner,
# options, gene list), copies this script there and submits the copy:
#   sbatch --output=OUT_DIR/slurm_logs/%A_%a.log [options] OUT_DIR/run/thoraxe_array.sbatch OUT_DIR
#
# Each task takes the genes TASK, TASK + COUNT, ... of genes.txt (largest exon table first, so
# the slow genes start first and land in different tasks) and runs, for each gene,
#   thoraxe -i IN_DIR/<gene> -o OUT_DIR/<gene> <options of run/thoraxe_args.txt>
# with the thoraxe package of run/code and the ProGraphMSA of run/bin. Each gene is logged in
# OUT_DIR/thoraxe_logs/<gene>.log, whose last lines are "resources={...}" (wall time, peak
# memory) and "exit=<code>". Genes with an exit line are skipped, so submitting again
# (thoraxe_array.py --resume OUT_DIR) resumes an interrupted run; a gene without one is
# rebuilt from scratch.
#
#SBATCH --job-name=thoraxe_array
#SBATCH --array=0-999
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=1-00:00:00
set -uo pipefail

if [ $# -ne 1 ]; then
  echo "usage: sbatch OUT_DIR/run/thoraxe_array.sbatch OUT_DIR (prepared by thoraxe_array.py)"
  exit 2
fi
OUT=$(readlink -f "$1")
RUN=$OUT/run
for f in input_dir python seed thoraxe_args.txt genes.txt code bin; do
  if [ ! -e "$RUN/$f" ]; then
    echo "missing $RUN/$f: prepare OUT_DIR with thoraxe_array.py"
    exit 1
  fi
done
IN=$(<"$RUN/input_dir")
PYTHON=$(<"$RUN/python")
SEED=$(<"$RUN/seed")
mapfile -t ARGS < "$RUN/thoraxe_args.txt"
mapfile -t GENES < <(grep -v '^[[:space:]]*$' "$RUN/genes.txt")
# an empty list is an error, not a task with nothing to do
if [ ${#GENES[@]} = 0 ]; then
  echo "no genes in $RUN/genes.txt"
  exit 1
fi

# a single BLAS thread: some numpy builds crash on import with threaded OpenBLAS
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# ThorAxe iterates over sets: the hash seed makes a rerun give the same output
export PYTHONHASHSEED=$SEED
export PATH=$RUN/bin:$PATH
unset PYTHONPATH
LOGS=$OUT/thoraxe_logs
mkdir -p "$LOGS" || exit 1

# runs the command given and prints its wall time and peak memory
TIMER='import json, resource, subprocess, sys, time
start = time.time()
code = subprocess.call(sys.argv[1:])
print("resources=" + json.dumps({
    "wall_s": round(time.time() - start, 1),
    "maxrss_mb": round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024, 1)}),
    flush=True)
sys.exit(code)'
# runs the thoraxe package of the snapshot given as first argument, and no other copy
# (first in sys.path: an editable install of thoraxe would win over PYTHONPATH)
THORAXE='import sys
code = sys.argv.pop(1)
sys.path.insert(0, code)
import thoraxe
if not thoraxe.__file__.startswith(code):
    sys.exit("thoraxe imported from %s, not from %s" % (thoraxe.__file__, code))
from thoraxe.thoraxe import main
sys.argv[0] = "thoraxe"
main()'

TASK=$(( ${SLURM_ARRAY_TASK_ID:-0} - ${SLURM_ARRAY_TASK_MIN:-0} ))
COUNT=${SLURM_ARRAY_TASK_COUNT:-1}
echo "task $TASK of $COUNT: genes $TASK, $((TASK + COUNT)), ... of ${#GENES[@]};" \
     "in=$IN host=$(hostname) start=$(date -Is)"

n_done=0 n_ok=0 n_skipped=0
for ((i = TASK; i < ${#GENES[@]}; i += COUNT)); do
  gene=${GENES[i]}
  log=$LOGS/$gene.log
  if grep -q '^exit=' "$log" 2>/dev/null; then
    n_skipped=$((n_skipped + 1))
    continue
  fi
  rm -rf "${OUT:?}/$gene"
  {
    echo "gene=$gene job=${SLURM_ARRAY_JOB_ID:-}_${SLURM_ARRAY_TASK_ID:-} host=$(hostname) start=$(date -Is)"
    "$PYTHON" -c "$TIMER" "$PYTHON" -c "$THORAXE" "$RUN/code/" \
        -i "$IN/$gene" -o "$OUT/$gene" "${ARGS[@]}"
  } > "$log" 2>&1
  rc=$?
  echo "exit=$rc" >> "$log"
  n_done=$((n_done + 1))
  [ "$rc" = 0 ] && n_ok=$((n_ok + 1))
done
echo "ran $n_done genes ($n_ok with exit 0), $n_skipped already done; end=$(date -Is)"
