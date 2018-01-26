#def get_deseq(wildcards):

def get_macs(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    sample_list = expand("macs/{sample}_peaks.narrowPeak", unit=unit.reset_index(), sample=wildcards.sample)
    return unit_list

rule merge_and_sort:
    input:
        "macs/{sample}_peaks.narrowPeak"
    output:
        "macs/bedgraph/{sample}_peaks.bedgraph"
    conda:
        "../envs/macs2.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} | bedtools merge -i - -c 5 -o mean > {output}
        """

rule unionbedgraph:
    input:
        expand("macs/bedgraph/{unit.sample}_peaks.bedgraph", unit=units.reset_index().drop_duplicates(subset='sample').itertuples())
    output:
        "macs/unionbedg/union.bedgraph"
    conda:
        "../envs/macs2.yaml"
    params:
        units=units
    shell:
        """
        bedtools unionbedg -i {input} > {output}.unsorted && sort -k1,1 -k2,2n {output}.unsorted | bedtools merge -i - > {output}
        """

rule multicov:
    input:
        bed="macs/unionbedg/union.bedgraph",
        bam=expand("dedup/{sample.sample}.bam", sample=samples.reset_index().itertuples())
    output:
        "counts/all_multicov.tsv"
    params:
        header="\t".join(samples.reset_index()['sample'].values)
    shell:
        """
        bedtools multicov -bams {input.bam} -bed {input.bed} > {output}
        """

rule coverage:
    input:
        bed="macs/unionbedg/union.bedgraph",
        bam="dedup/{sample}-{unit}.bam"
    conda:
        "../envs/macs2.yaml"
    output:
        "macs/coverage/{sample}-{unit}-peaks_coverage.bedgraph"
    params: ""
    shell:
        """
        bedtools multicov -bams {input.bam} -bed {input.bed} > {output}
        """

rule count_matrix:
    input:
        expand("macs/coverage/{sample.sample}-peaks_coverage.bedgraph", sample=samples.reset_index().itertuples())
    output:
        "counts/all.tsv"
    params:
        samples=samples
    script:
        "../scripts/count-matrix.py"

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        counts="counts/all.tsv"
    output:
        "deseq2/all.rds"
    params:
        samples=config["samples"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads(None)
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "deseq2/all.rds"
    output:
        "results/pca.svg"
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]

rule deseq2_plot:
    input:
        "deseq2/all.rds"
    output:
        volcano="results/diffexp/{contrast}.volcano.pdf",
        ma_plot="results/diffexp/{contrast}.ma-plot.pdf",
    params:
        contrast=get_contrast,
        labels=config["diffexp"]["label"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"

rule deseq2_write:
    input:
        "deseq2/all.rds"
    output:
        table="results/diffexp/{contrast}.diffexp.tsv",
        sig="results/diffexp/{contrast}_sig.diffexp.tsv",
        gain="results/diffexp/{contrast}_gain.diffexp.tsv",
        loss="results/diffexp/{contrast}_loss.diffexp.tsv",
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.write.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2-write.R"

rule deseq2_to_bed:
    input:
        gain="results/diffexp/{contrast}_gain.diffexp.tsv",
        loss="results/diffexp/{contrast}_loss.diffexp.tsv",
    output:
        gain_bed="results/diffexp/bed/{contrast}_gain.diffexp.bed",
        loss_bed="results/diffexp/bed/{contrast}_loss.diffexp.bed",
    shell:
        """
        tail -n +2 {input.gain} | cut -d' ' -f1 | cut -d'_' -f1,2,3 --output-delimiter=$'\t' > {output.gain_bed} && tail -n +2 {input.loss} | cut -d' ' -f1 | cut -d'_' -f1,2,3 --output-delimiter=$'\t' > {output.loss_bed}
        """

rule homer:
    input:
        gain="results/diffexp/bed/{contrast}_gain.diffexp.bed",
        loss="results/diffexp/bed/{contrast}_loss.diffexp.bed",
    output:
        gain="motif/{contrast}_gain_sites/homerResults.html",
        loss="motif/{contrast}_loss_sites/homerResults.html",
    conda:
        "../envs/homer.yaml"
    params:
        gain="motif/{contrast}_gain_sites/",
        loss="motif/{contrast}_loss_sites/",
        genome="/mnt/data1/John/bin/pyenv/versions/anaconda3-5.0.1/pkgs/homer-4.9.1-pl5.22.0_5/share/homer-4.9.1-5/data/genomes/mm10"
    shell:
        """
        findMotifsGenome.pl {input.gain} {params.genome} {params.gain} -p 30 && findMotifsGenome.pl {input.loss} {params.genome} {params.loss} -p 30
        """
