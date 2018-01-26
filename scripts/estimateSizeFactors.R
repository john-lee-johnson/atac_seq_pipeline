log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

coverage <- read.table(snakemake@input[[1]], header=TRUE, row.names="peak", check.names = FALSE)
sizeFactors <- estimateSizeFactorsForMatrix(coverage)

write.table(sizeFactors, file=snakemake@output[[1]], quote = FALSE, col.names = FALSE)
