#!/bin/bash
# Download and index the Ensembl bulk files of one release with ensembl_bulk, on SLURM.
#
#   sbatch [sbatch options] scripts/ensembl_bulk.sbatch RELEASE PREFIX [download options]
#
# The files are downloaded into PREFIX_downloading, which is renamed PREFIX_<YYYYMMDD>
# (the day the download finishes) before the indexes are built. Submitting the same
# command again resumes: verified files and built indexes are kept, and an already
# renamed folder is indexed directly. Extra arguments go to "ensembl_bulk download",
# e.g. "-l species.txt" to download only some species (default: every species with an
# ortholog of the query species; about 120 GB for human in Ensembl 115-116, and about
# 0.55 TB once indexed).
#
# Environment: PYTHON, the Python interpreter where thoraxe is installed (default: python).
#
# Example:
#   sbatch --output=ensembl116_%j.log scripts/ensembl_bulk.sbatch 116 /data/ensembl116
#
#SBATCH --job-name=ensembl_bulk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
set -uo pipefail

if [ $# -lt 2 ]; then
  echo "usage: sbatch scripts/ensembl_bulk.sbatch RELEASE PREFIX [download options]"
  exit 2
fi
RELEASE=$1
PREFIX=$2
shift 2
PYTHON=${PYTHON:-python}
# a single BLAS thread: some numpy builds crash on import with threaded OpenBLAS
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
BULK=("$PYTHON" -m thoraxe.bulk_ensembl.cli)

echo "Ensembl $RELEASE -> ${PREFIX}_<date>; $("${BULK[@]}" --version); host=$(hostname) start=$(date -Is)"

FINAL=$(ls -d "${PREFIX}"_2[0-9][0-9][0-9][0-9][0-9][0-9][0-9] 2>/dev/null | tail -n 1)
if [ -z "$FINAL" ]; then
  "${BULK[@]}" download "${PREFIX}_downloading" --release "$RELEASE" "$@" \
    || { echo "download FAILED"; exit 1; }
  FINAL=${PREFIX}_$(date +%Y%m%d)
  mv "${PREFIX}_downloading" "$FINAL" || exit 1
  echo "download done $(date -Is) -> $FINAL"
fi
"${BULK[@]}" index "$FINAL" --workers "${SLURM_CPUS_PER_TASK:-8}" || { echo "index FAILED"; exit 1; }
du -sh "$FINAL"
echo "end=$(date -Is)"
