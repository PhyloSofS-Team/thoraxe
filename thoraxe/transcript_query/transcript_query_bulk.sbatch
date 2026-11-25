#!/bin/bash
# Build the transcript_query folders of many genes from an Ensembl bulk folder, on SLURM.
#
#   sbatch [sbatch options] thoraxe/transcript_query/transcript_query_bulk.sbatch \
#          BULK_DIR OUT_DIR [GENE_LIST] [transcript_query options]
#
# Writes OUT_DIR/<gene>/Ensembl/ for every gene of GENE_LIST (one Ensembl gene id per line),
# with "transcript_query <gene> --bulk-dir BULK_DIR --orthology 1:n" followed by the
# transcript_query options given (e.g. "-l species.txt"). Without GENE_LIST, the genes are
# every protein-coding gene of the query species of BULK_DIR, listed in OUT_DIR/gene_list.txt.
# BULK_DIR is a folder prepared with ensembl_bulk (thoraxe/bulk_ensembl/ensembl_bulk.sbatch).
#
# The genes are spread over the array tasks (100 by default, change it with --array). A gene
# takes about 10 s with 19 species and about 1 minute with all ~200 species (peak memory
# 1.5 GB), so with 100 tasks the 20,000 human protein-coding genes take about 30 minutes or
# about 4 hours per task; use more tasks (e.g. --array=0-399) to finish sooner. Each gene is logged in
# OUT_DIR/transcript_query_logs/<gene>.log, whose last line is "exit=<code>". Genes with
# such a line are skipped, so submitting the same command again resumes an interrupted run;
# a gene without it is rebuilt from scratch. Genes without orthologs end with exit=1, as
# with the online transcript_query, and keep only Ensembl/ensembl_version.csv (and tree.nh).
#   failed or without orthologs:  grep -L '^exit=0' OUT_DIR/transcript_query_logs/*.log
#
# Environment: PYTHON, the Python interpreter where thoraxe is installed (default: python).
#
# Example:
#   sbatch --output=genes_%A_%a.log thoraxe/transcript_query/transcript_query_bulk.sbatch \
#          /data/ensembl116_20261125 /data/ensembl116_genes
#
#SBATCH --job-name=transcript_query_bulk
#SBATCH --array=0-99
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1-00:00:00
set -uo pipefail

if [ $# -lt 2 ]; then
  echo "usage: sbatch thoraxe/transcript_query/transcript_query_bulk.sbatch BULK_DIR OUT_DIR [GENE_LIST] [transcript_query options]"
  exit 2
fi
BULK=$(readlink -f "$1")
OUT=$2
shift 2
GENE_LIST=
if [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; then
  GENE_LIST=$(readlink -f "$1")
  shift
fi
PYTHON=${PYTHON:-python}
# a single BLAS thread: some numpy builds crash on import with threaded OpenBLAS
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}

mkdir -p "$OUT/transcript_query_logs" || exit 1
OUT=$(readlink -f "$OUT")
LOGS=$OUT/transcript_query_logs

if [ -z "$GENE_LIST" ]; then
  # every protein-coding gene of the query species, in a stable order (same list in every task)
  GENE_LIST=$OUT/gene_list.txt
  "$PYTHON" - "$BULK" > "$GENE_LIST.$$" <<'EOF' || { echo "cannot list the genes of $BULK"; exit 1; }
import sqlite3
import sys
from thoraxe.bulk_ensembl.layout import BulkDir
bulk = BulkDir(sys.argv[1])
con = sqlite3.connect(bulk.gtf_db(bulk.query_species))
for (gene,) in con.execute(
        "SELECT gene_id FROM gene WHERE biotype = 'protein_coding' ORDER BY gene_id"):
    print(gene)
EOF
  mv "$GENE_LIST.$$" "$GENE_LIST"
fi
mapfile -t GENES < <(grep -v '^[[:space:]]*$' "$GENE_LIST")

TASK=$(( ${SLURM_ARRAY_TASK_ID:-0} - ${SLURM_ARRAY_TASK_MIN:-0} ))
COUNT=${SLURM_ARRAY_TASK_COUNT:-1}
echo "task $TASK of $COUNT: genes $TASK, $((TASK + COUNT)), ... of ${#GENES[@]} from $GENE_LIST;" \
     "bulk=$BULK host=$(hostname) start=$(date -Is)"

n_done=0 n_ok=0 n_skipped=0
for ((i = TASK; i < ${#GENES[@]}; i += COUNT)); do
  gene=${GENES[i]}
  log=$LOGS/$gene.log
  if grep -q '^exit=' "$log" 2>/dev/null; then
    n_skipped=$((n_skipped + 1))
    continue
  fi
  rm -rf "${OUT:?}/$gene"
  (cd "$OUT" && "$PYTHON" -m thoraxe.transcript_query.transcript_query "$gene" \
      --bulk-dir "$BULK" --orthology 1:n "$@") > "$log" 2>&1
  rc=$?
  echo "exit=$rc" >> "$log"
  n_done=$((n_done + 1))
  [ "$rc" = 0 ] && n_ok=$((n_ok + 1))
done
echo "built $n_done genes ($n_ok with exit 0), $n_skipped already done; end=$(date -Is)"
