#!/bin/bash
#SBATCH --job-name=fastqc_multiqc
#SBATCH --cpus-per-task=12
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --output=/home/rona/mvirag/Paryani_2024/Analysis/Outputs/qc_%j.out

set -euo pipefail

export PATH="$HOME/software/FastQC:$HOME/software/multiqc_env/bin:$PATH"

BASE_DIR="/home/rona/mvirag/Paryani_2024"
QCSET="${BASE_DIR}/Data/QC_ready"
FASTQC_OUT="${BASE_DIR}/Analysis/qc_reports/fastqc"
MULTIQC_OUT="${BASE_DIR}/Analysis/qc_reports/multiqc"

command -v fastqc  >/dev/null 2>&1 || { echo "HIBA: fastqc nem található."; exit 1; }
command -v multiqc >/dev/null 2>&1 || { echo "HIBA: multiqc nem található."; exit 1; }

mkdir -p "$FASTQC_OUT" "$MULTIQC_OUT"

N=$(ls "${QCSET}"/*.fastq.gz | wc -l)
echo "--- FastQC indul ${N} fájlon ---"

fastqc -t 6 -o "$FASTQC_OUT" "${QCSET}"/*.fastq.gz

echo "--- MultiQC összesítés ---"
multiqc "$FASTQC_OUT" -o "$MULTIQC_OUT" --force

echo "Kész. HTML: ${MULTIQC_OUT}/multiqc_report.html"

