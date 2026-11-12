#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --array=1-325%8
#SBATCH --output=/home/rona/mvirag/Paryani_2024/Analysis/Outputs/sra_%A_%a.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=10:00:00

set -euo pipefail

# ============================================================
# EZEKET KELL ÁTÍRNI ÚJ CIKKNÉL / PROJEKTNÉL:
# ============================================================
BASE_DIR="/home/rona/mvirag/Paryani_2024"
ACC_LIST="${BASE_DIR}/Data/Metadata/SRR_Acc_List.txt"
# ============================================================

export PATH="/home/rona/mvirag/software/sratoolkit.3.4.1-ubuntu64/bin:$PATH"

command -v prefetch >/dev/null 2>&1 || { echo "HIBA: prefetch nem található."; exit 1; }
command -v fasterq-dump >/dev/null 2>&1 || { echo "HIBA: fasterq-dump nem található."; exit 1; }

OUTDIR="${BASE_DIR}/Data/Raw_fastq"
DONE_MARKER_DIR="${BASE_DIR}/Data/Metadata/.done_markers"

mkdir -p "$OUTDIR" "$DONE_MARKER_DIR"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACC_LIST" | tr -d '\r' | xargs)
[[ -z "$SRR" ]] && { echo "Üres sor a(z) ${SLURM_ARRAY_TASK_ID}. indexen."; exit 0; }

DONE_MARKER="${DONE_MARKER_DIR}/${SRR}.done"
if [[ -f "$DONE_MARKER" ]]; then
    echo "Már kész, kihagyva: $SRR"
    exit 0
fi

echo "--- Feldolgozás: $SRR (array index: ${SLURM_ARRAY_TASK_ID}) ---"

DF_AVAIL=$(df --output=avail "$OUTDIR" | tail -1)
if [[ "$DF_AVAIL" -lt 20000000 ]]; then
    echo "HIBA: kevesebb, mint 20GB szabad hely — leállítás."
    exit 1
fi

MAX_RETRIES=3
for attempt in $(seq 1 $MAX_RETRIES); do
    if prefetch "$SRR" -O "$OUTDIR"; then
        break
    fi
    echo "Próba $attempt sikertelen: $SRR"
    if [[ "$attempt" -eq "$MAX_RETRIES" ]]; then
        echo "HIBA: végleges sikertelenség $MAX_RETRIES próba után: $SRR"
        exit 1
    fi
    sleep $((attempt * 30))
done

fasterq-dump --split-files --outdir "$OUTDIR" --threads 2 "${OUTDIR}/${SRR}/${SRR}.sra"

shopt -s nullglob
for f in "${OUTDIR}/${SRR}"*.fastq; do
    gzip -f "$f"
done
shopt -u nullglob

mkdir -p "${BASE_DIR}/Data/Stale_sra"
mv "${OUTDIR:?}/${SRR:?}" "${BASE_DIR}/Data/Stale_sra/"

touch "$DONE_MARKER"
echo "Sikeresen kész: $SRR"
