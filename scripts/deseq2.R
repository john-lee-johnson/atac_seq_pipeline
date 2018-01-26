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

# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, main=snakemake@params[["labels"]], ylim=c(-2,2))
dev.off()

resSig.gain = res[ res$padj < 0.05 & res$log2FoldChange > 0.55 , ]
resSig.loss = res[ res$padj < 0.05 & res$log2FoldChange < 0.55 , ]
resSig = res[ res$padj < 0.05 & abs(res$log2FoldChange) > 0.55 , ]

pdf(snakemake@output[["volcano"]])
#     width = 480, height = 480, units = "px", pointsize = 12,
#     compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
#     bg = "white", res = NA,  ...,
#     type = c("cairo", "Xlib", "quartz"), antialias)

# Make a basic volcano plot

with(as.data.frame(res), plot(log2FoldChange, -log10(pvalue), pch=20, main=snakemake@params[["labels"]], xlim=c(-2.5,2.5), , ylim=c(0,40)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(as.data.frame(res), padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(as.data.frame(res), abs(log2FoldChange)>0.55), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(as.data.frame(res), padj<.05 & abs(log2FoldChange)>0.55), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
dev.off()
