#!/bin/bash
set -uo pipefail

RAW="/home/rona/mvirag/Paryani_2024/Data/Raw_fastq"
META="/home/rona/mvirag/Paryani_2024/Data/Metadata"
QCSET="/home/rona/mvirag/Paryani_2024/Data/QC_ready"

rm -rf "$QCSET"
mkdir -p "$QCSET"
> "${META}/qc_set_manifest.txt"
> "${META}/qc_set_duplicates.txt"

while read -r SRR; do
    SRR=$(echo "$SRR" | tr -d '\r' | xargs)
    [[ -z "$SRR" ]] && continue

    shopt -s nullglob
    FILES=("${RAW}/${SRR}"*.fastq.gz)
    shopt -u nullglob
    [[ ${#FILES[@]} -eq 0 ]] && continue

    NCBI=()
    ENA=()
    for f in "${FILES[@]}"; do
        HDR=$(zcat "$f" 2>/dev/null | head -1)
        if [[ "$HDR" == *"length="* ]]; then
            NCBI+=("$f")
        elif [[ "$HDR" =~ /[0-9]$ ]]; then
            ENA+=("$f")
        else
            NCBI+=("$f")
        fi
    done

    if [[ ${#NCBI[@]} -gt 0 && ${#ENA[@]} -gt 0 ]]; then
        CHOSEN=("${NCBI[@]}")
        echo "${SRR}: NCBI (${#NCBI[@]} fájl) megtartva, ENA (${#ENA[@]} fájl) kihagyva" >> "${META}/qc_set_duplicates.txt"
    elif [[ ${#NCBI[@]} -gt 0 ]]; then
        CHOSEN=("${NCBI[@]}")
    else
        CHOSEN=("${ENA[@]}")
    fi

    for f in "${CHOSEN[@]}"; do
        ln -s "$f" "${QCSET}/$(basename "$f")"
        echo "$(basename "$f")" >> "${META}/qc_set_manifest.txt"
    done

done < "${META}/list_complete.txt"

echo "=== QC KÉSZLET ==="
echo "Minta:               $(wc -l < ${META}/list_complete.txt)"
echo "Fájl a QC_ready-ben: $(ls ${QCSET} | wc -l)"
echo "Vegyes forrású minta: $(wc -l < ${META}/qc_set_duplicates.txt)"
echo "Törött link:         $(find -L ${QCSET} -type l | wc -l)"
