#!/bin/bash
#SBATCH --job-name=salmon_quant_tseng
#SBATCH --output=/home/rona/mvirag/arcticles/Tseng_2016/logs/salmon_%A_%a.log
#SBATCH --array=1-14
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --partition=rall
#SBATCH --nodelist=quad-5,quad-6,twin-1,twin-2,twin-3,twin-4

set -euo pipefail
export PATH="$HOME/.cargo/bin:$PATH"

BASE="/home/rona/mvirag/arcticles/Tseng_2016"
TRIMDIR="$BASE/Data/Trimmed_fastq"
OUTDIR="$BASE/Analysis/Outputs/salmon_quant"
INDEX="/home/rona/mvirag/arcticles/Reference/gencode_v44_salmon_index"
ACCLIST="$BASE/Data/Metadata/SRR_Acc_List.txt"

mkdir -p "$OUTDIR"

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCLIST")
R1="${TRIMDIR}/${SRR}_R1.trim.fastq.gz"
R2="${TRIMDIR}/${SRR}_R2.trim.fastq.gz"
SINGLE="${TRIMDIR}/${SRR}.trim.fastq.gz"

if [[ -f "$R1" && -f "$R2" ]]; then
    echo "Futtatás [Paired-End]: $SRR"
    salmon quant -i "$INDEX" -l A -1 "$R1" -2 "$R2" -p 8 \
                 --validateMappings -o "${OUTDIR}/${SRR}_quant"
elif [[ -f "$SINGLE" ]]; then
    echo "Futtatás [Single-End]: $SRR"
    salmon quant -i "$INDEX" -l A -r "$SINGLE" -p 8 \
                 --validateMappings -o "${OUTDIR}/${SRR}_quant"
else
    echo "Hiba: A ${SRR} mintához tartozó trimmelt fájl(ok) nem találhatóak!"
    exit 1
fi
