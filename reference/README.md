# Reference files (GENCODE release 44)

Both datasets were quantified against the same Salmon index and summarised to genes with
the same annotation release.

## Salmon index
- Built from `gencode.v44.transcripts.fa.gz` with Salmon 2.4.1: k-mer size 31, no decoy
  sequences, 251,955 reference sequences.
- The build command was not saved; the settings are taken from the index's own record,
  stored here as `salmon_index_info.json`.

## tx2gene_v44.csv
- 252,834 transcripts with their gene (columns `transcript_id`, `gene_id`), extracted from
  `gencode.v44.annotation.gtf` by `make_tx2gene_v44.R`.
- Used for Paryani et al. (2024). For Lin et al. (2016) the same transcript-to-gene links were
  read from the transcript names in the Salmon output (see `scripts/Lin_2016/05_tximport_counts.R`).

## Source files (not stored here)
Download from GENCODE release 44: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/
- `gencode.v44.transcripts.fa.gz`
- `gencode.v44.annotation.gtf.gz`
