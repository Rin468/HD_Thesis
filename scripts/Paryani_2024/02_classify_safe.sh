#!/bin/bash
set -uo pipefail

RAW="/home/rona/mvirag/Paryani_2024/Data/Raw_fastq"
META="/home/rona/mvirag/Paryani_2024/Data/Metadata"

> "${META}/list_complete.txt"
> "${META}/list_incomplete.txt"
> "${META}/files_to_reject.txt"
> "${META}/list_skipped.txt"

while read -r SRR; do
    SRR=$(echo "$SRR" | tr -d '\r' | xargs)
    [[ -z "$SRR" ]] && continue

    # 1) kihagyjuk, ha a job éppen dolgozik rajta
    if grep -qx "$SRR" "${META}/in_progress.txt" 2>/dev/null; then
        echo "${SRR}  (éppen fut)" >> "${META}/list_skipped.txt"
        continue
    fi

    # 2) kihagyjuk, ha van aktív .sra.lock (letöltés alatt)
    if [[ -e "${RAW}/${SRR}/${SRR}.sra.lock" ]]; then
        echo "${SRR}  (sra.lock aktív)" >> "${META}/list_skipped.txt"
        continue
    fi

    shopt -s nullglob
    FILES=("${RAW}/${SRR}"*.fastq.gz)
    shopt -u nullglob

    [[ ${#FILES[@]} -eq 0 ]] && continue

    BAD=0
    for f in "${FILES[@]}"; do
        N=$(zcat "$f" 2>/dev/null | head -400 | awk 'NR%4==0' | fold -w1 | sort -u | wc -l)
        if [[ "$N" -le 1 ]]; then
            echo "$f" >> "${META}/files_to_reject.txt"
            BAD=$((BAD+1))
        fi
    done

    if [[ "$BAD" -eq 0 ]]; then
        echo "$SRR" >> "${META}/list_complete.txt"
    else
        echo "$SRR" >> "${META}/list_incomplete.txt"
    fi

done < "${META}/SRR_Acc_List.txt"

echo "=== EREDMÉNY ==="
echo "Hibátlan, QC-re kész:  $(wc -l < ${META}/list_complete.txt)"
echo "Részben hibás:         $(wc -l < ${META}/list_incomplete.txt)"
echo "Kihagyva (fut):        $(wc -l < ${META}/list_skipped.txt)"
echo "Áthelyezendő fájl:     $(wc -l < ${META}/files_to_reject.txt)"
