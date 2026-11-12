#!/bin/bash
#SBATCH --job-name=salmon_quant_smart
#SBATCH --output=/home/rona/mvirag/arcticles/Paryani_2024/logs/salmon_smart_%A_%a.log
#SBATCH --array=1-325%32
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --partition=bquads
#SBATCH --nodelist=quad-6

set -euo pipefail

# Salmon elérési út beállítása
export PATH="$HOME/.cargo/bin:$PATH"

BASE="/home/rona/mvirag/arcticles/Paryani_2024"
TRIMDIR="$BASE/Data/Trimmed_fastq"
OUTDIR="$BASE/Analysis/Outputs/salmon_quant"
INDEX="/home/rona/mvirag/arcticles/Reference/gencode_v44_salmon_index"
ACCLIST="$BASE/Data/Metadata/SRR_Acc_List.txt"

mkdir -p "$OUTDIR"

# Az aktuális SRR azonosító lekérése a tömbből
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCLIST")
R1="${TRIMDIR}/${SRR}_R1.trim.fastq.gz"
R2="${TRIMDIR}/${SRR}_R2.trim.fastq.gz"
SINGLE="${TRIMDIR}/${SRR}.trim.fastq.gz"

# 1. eset: Paired-end fájlok léteznek
if [[ -f "$R1" && -f "$R2" ]]; then
    echo "Futtatás [Paired-End]: $SRR"
    salmon quant -i "$INDEX" \
                 -l A \
                 -1 "$R1" \
                 -2 "$R2" \
                 -p 8 \
                 --validateMappings \
                 -o "${OUTDIR}/${SRR}_quant"

# 2. eset: Single-end fájl létezik
elif [[ -f "$SINGLE" ]]; then
    echo "Futtatás [Single-End]: $SRR"
    salmon quant -i "$INDEX" \
                 -l A \
                 -r "$SINGLE" \
                 -p 8 \
                 --validateMappings \
                 -o "${OUTDIR}/${SRR}_quant"

# 3. eset: Egyik sem található
else
    echo "Hiba: A ${SRR} mintához tartozó trimmelt fájl(ok) nem találhatóak!"
    exit 1
fi
