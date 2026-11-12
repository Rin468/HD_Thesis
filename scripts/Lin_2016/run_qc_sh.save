#!/bin/bash
#SBATCH --job-name=tseng_qc
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/home/rona/mvirag/arcticles/Tseng_2016/Analysis/Outputs/qc_%j.out

set -euo pipefail

# Szoftver útvonalak betöltése a szerver környezethez
export PATH="$HOME/software/FastQC:$HOME/software/multiqc_env/bin:$PATH"

BASE_DIR="/home/rona/mvirag/arcticles/Tseng_2016"
FASTQ_DIR="${BASE_DIR}/Data/Raw_fastq"
FASTQC_OUT="${BASE_DIR}/Analysis/qc_reports/fastqc"
MULTIQC_OUT="${BASE_DIR}/Analysis/qc_reports/multiqc"

# Ellenőrizzük, hogy a programok elérhetők-e
command -v fastqc  >/dev/null 2>&1 || { echo "HIBA: fastqc nem található."; exit 1; }
command -v multiqc >/dev/null 2>&1 || { echo "HIBA: multiqc nem található."; exit 1; }

# Mappák létrehozása, ha nem léteznek
mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

N=$(ls "${FASTQ_DIR}"/*.fastq.gz | wc -l)
echo "--- FastQC indul ${N} fájlon ---"

# FastQC futtatása párhuzamos szállakkal
fastqc -t 6 -o "$FASTQC_OUT" "${FASTQ_DIR}"/*.fastq.gz

echo "--- MultiQC összesítés ---"
multiqc "$FASTQC_OUT" -o "$MULTIQC_OUT" --force

echo "Kész. HTML: ${MULTIQC_OUT}/multiqc_report.html"
