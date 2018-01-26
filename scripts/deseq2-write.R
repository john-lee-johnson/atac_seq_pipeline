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

dds <- readRDS(snakemake@input[[1]])
contrast <- c("condition", snakemake@params[["contrast"]])

res <- results(dds, contrast=contrast, parallel=parallel)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)
# sort by p-value
res <- res[order(res$padj),]
res <- res[order(-res$log2FoldChange),]
# TODO explore IHW usage

# Remove NA
res <- res[!is.na(res$padj), ]
res <- res[!is.na(res$log2FoldChange), ]

# Threshold based on log2FoldChange and padj
resSig.gain = res[ res$padj < 0.05 & res$log2FoldChange > 0.55 , ]
resSig.loss = res[ res$padj < 0.05 & res$log2FoldChange < 0.55 , ]
resSig = res[ res$padj < 0.05 & abs(res$log2FoldChange) > 0.55 , ]

# Write out to file to convert to bed
write.table(as.data.frame(resSig.gain), file=snakemake@output[["gain"]], quote = FALSE)
write.table(as.data.frame(resSig.loss), file=snakemake@output[["loss"]], quote = FALSE)
write.table(as.data.frame(resSig), file=snakemake@output[["sig"]], quote = FALSE)
write.table(as.data.frame(res), file=snakemake@output[["table"]], quote = FALSE)
