#!/bin/bash
#SBATCH --job-name=trim_fastp_fast
#SBATCH --array=1-325%24
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=03:00:00
#SBATCH --output=/home/rona/mvirag/arcticles/Paryani_2024/logs/trim_%A_%a.log
#SBATCH --partition=bquads
#SBATCH --nodelist=quad-6

set -euo pipefail

fastp="/home/rona/mvirag/tools/fastp"

BASE="/home/rona/mvirag/arcticles/Paryani_2024"
RAWDIR="$BASE/Data/Raw_fastq"
TRIMDIR="$BASE/Data/Trimmed_fastq"
REPDIR="$BASE/Analysis/qc_reports/fastp"
DONEDIR="$BASE/Data/Metadata/.done_markers/trim"
ACCLIST="$BASE/Data/Metadata/SRR_Acc_List.txt"

mkdir -p "$TRIMDIR" "$REPDIR" "$DONEDIR"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCLIST")
DONEFLAG="${DONEDIR:?}/${SRR:?}.done"
[[ -f "$DONEFLAG" ]] && { echo "$SRR already trimmed, skipping"; exit 0; }

# generikus párkeresés: az összes SRR-hez tartozó fastq.gz fájl begyűjtése,
# numerikus suffix szerint rendezve
mapfile -t FILES < <(ls "${RAWDIR}/${SRR}"*.fastq.gz 2>/dev/null | sort -V)

case "${#FILES[@]}" in
    2)
        R1="${FILES[0]}"
        R2="${FILES[1]}"
        "$fastp" \
            -i "$R1" -I "$R2" \
            -o "${TRIMDIR}/${SRR}_R1.trim.fastq.gz" -O "${TRIMDIR}/${SRR}_R2.trim.fastq.gz" \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 25 \
            --thread 8 \
            --json "${REPDIR}/${SRR}.fastp.json" \
            --html "${REPDIR}/${SRR}.fastp.html"
        ;;
    1)
        SE="${FILES[0]}"
        "$fastp" \
            -i "$SE" \
            -o "${TRIMDIR}/${SRR}.trim.fastq.gz" \
            --qualified_quality_phred 20 \
            --length_required 25 \
            --thread 8 \
            --json "${REPDIR}/${SRR}.fastp.json" \
            --html "${REPDIR}/${SRR}.fastp.html"
        ;;
    *)
        echo "ERROR: ${#FILES[@]} fastq.gz file(s) found for $SRR in $RAWDIR (expected 1 or 2): ${FILES[*]}" >&2
        exit 1
        ;;
esac

touch "$DONEFLAG"
