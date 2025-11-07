#!/usr/bin/env Rscript
# Peak annotation using ChIPseeker

suppressMessages({
    library(ChIPseeker)
    library(GenomicFeatures)
})

# Get input arguments from Snakemake
peak_file <- snakemake@input$peaks
gtf_file <- snakemake@input$gtf
output_file <- snakemake@output[[1]]

# Determine genome build from config or filename
genome_build <- ifelse(grepl("mm10|mouse", basename(gtf_file), ignore.case=TRUE), "mm10", "hg38")

# Load appropriate TxDb
if (genome_build == "mm10") {
    suppressMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
    suppressMessages(library(org.Mm.eg.db))
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    orgdb <- org.Mm.eg.db
} else {
    suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
    suppressMessages(library(org.Hs.eg.db))
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    orgdb <- org.Hs.eg.db
}

cat("Annotating peaks for", genome_build, "genome...\n")

# Read peaks
peaks <- readPeakFile(peak_file)

# Annotate peaks
peak_anno <- annotatePeak(peaks, 
                         tssRegion=c(-3000, 3000),
                         TxDb=txdb, 
                         annoDb=orgdb)

# Convert to data frame and save
peak_anno_df <- as.data.frame(peak_anno)
write.table(peak_anno_df, output_file, sep="\t", row.names=FALSE, quote=FALSE)

cat("Peak annotation completed. Results saved to:", output_file, "\n")
cat("Total peaks annotated:", nrow(peak_anno_df), "\n")

# Print summary
cat("\nAnnotation summary:\n")
print(table(peak_anno_df$annotation))