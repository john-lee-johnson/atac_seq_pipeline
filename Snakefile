import pandas as pd
shell.executable("bash")

configfile: "config.yaml"
samples = pd.read_table(config["samples"], index_col="sample")
units = pd.read_table(config["units"], index_col=["sample", "lane"], dtype=str).sort_index(level=0)

wildcard_constraints:
    sample="\w+",
    unit="\w+",
    lane="\w+",
    fastq="\w+"

def is_single_end(sample, lane):
    return pd.isnull(units.loc[(sample, lane), "fq2"])

rule all:
    input:
        expand("idr/{unit.sample}/{unit.sample}-optimal-peakset.narrowPeak", unit=units.reset_index().drop_duplicates(subset='sample').itertuples()),
        expand("idr/{unit.sample}/{unit.sample}-conservative-peakset.narrowPeak", unit=units[units['unit'].notnull()].reset_index().itertuples()),
        expand("macs/coverage/{unit.sample}-{unit.unit}-peaks_coverage.bedgraph", unit=units[units['unit'].notnull()].reset_index().itertuples()),
        expand("macs/coverage/{unit.sample}-peaks_coverage.bedgraph", unit=units.reset_index().drop_duplicates(subset='sample').itertuples()),
        expand("bw/rpkm/units/{unit.sample}-{unit.unit}-RPKM.bw", unit=units[units['unit'].notnull()].reset_index().itertuples()),
        expand("bw/rpkm/{unit.sample}-RPKM.bw", unit=units.reset_index().drop_duplicates(subset='sample').itertuples()),
        expand("bw/macs2/units/{unit.sample}-{unit.unit}_treat_pileup.bw", unit=units[units['unit'].notnull()].reset_index().itertuples()),
        expand("bw/macs2/{unit.sample}_treat_pileup.bw", unit=units.reset_index().drop_duplicates(subset='sample').itertuples()),
        expand("results/diffexp/{contrast}.diffexp.tsv",
            contrast=config["diffexp"]["contrasts"]),
        expand("results/pca.svg",
            contrast=config["diffexp"]["contrasts"]),
        expand("results/diffexp/{contrast}.volcano.pdf",
            contrast=config["diffexp"]["contrasts"]),
        expand("motif/{contrast}_{results}_sites/homerResults.html",
            contrast=config["diffexp"]["contrasts"], results=['gain', 'loss'])

rule qc:
    input:
        expand("qc/{fastq}_fastqc.html", fastq=units[['fq1', 'fq2']].stack().str.replace('.fastq.gz','')),
        expand("qc/{unit.sample}-{unit.unit}-{unit.lane}_{group}_trimmed_fastqc.html", unit=units[units.notnull()['unit']].reset_index().set_index(['sample', 'unit', 'lane']).sort_index(level=0).stack().reset_index().itertuples(), group=[1,2]),
        expand("qc/{unit.sample}-{unit.lane}_{group}_trimmed_fastqc.html", unit=units[units.isnull()['unit']].reset_index().set_index(['sample', 'lane']).sort_index(level=0).stack().reset_index().itertuples(), group=[1,2]),
        expand("stats/{unit.sample}-{unit.unit}.isize.pdf", unit=units[units['unit'].notnull()].reset_index().itertuples()),
        expand("stats/{unit.sample}.isize.pdf", unit=units.reset_index().drop_duplicates(subset='sample').itertuples()),

include: "rules/fastqc.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/idr.smk"
include: "rules/macs2.smk"
include: "rules/merge.smk"
include: "rules/blacklist.smk"
include: "rules/mark_duplicates.smk"
include: "rules/bamtobedpe.smk"
include: "rules/pseudoreplicates.smk"
include: "rules/make_bigwig.smk"
include: "rules/coverage.smk"
include: "rules/homer.smk"
include: "rules/insert_size.smk"
include: "rules/diffexp.smk"
