# Build the transcript-to-gene table for GENCODE release 44.
# Same code as in Paryani_DE_analysis.R (server), where tx2gene_v44.csv was created;
# only the file paths were simplified.
# Input: gencode.v44.annotation.gtf (GENCODE release 44, unzipped)

library(rtracklayer)

gtf <- import("gencode.v44.annotation.gtf")
tx2gene <- as.data.frame(
  mcols(gtf)[gtf$type == "transcript", c("transcript_id", "gene_id")]
)
write.csv(tx2gene, "tx2gene_v44.csv", row.names = FALSE)
