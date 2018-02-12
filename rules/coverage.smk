def get_bam_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
    # has units
        return "bam/merge/{sample}.bam".format(sample=wildcards.sample)
    else:
    # no units
        return "bam/{sample}.bam".format(sample=wildcards.sample)

def get_idr_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
    # has units
        return "idr/{sample}/{sample}-conservative-peakset.narrowPeak".format(sample=wildcards.sample)
    else:
    # no units
        return "idr/{sample}/{sample}-optimal-peakset.narrowPeak".format(sample=wildcards.sample)

def get_idr_subset_input(wildcards):
    units = pd.read_table(config["subset"], index_col=["sample"], dtype=str).sort_index(level=0)
    idr_list = []
    idr_list = idr_list + expand("idr/{unit.sample}/{unit.sample}-conservative-peakset.narrowPeak", unit=units[units['unit'].notnull()].reset_index().drop_duplicates(subset='sample').itertuples())
    idr_list = idr_list + expand("idr/{unit.sample}/{unit.sample}-optimal-peakset.narrowPeak", unit=units[units['unit'].isnull()].reset_index().drop_duplicates(subset='sample').itertuples())
    return(idr_list)

def get_idr_subset_input_distal(wildcards):
    units = pd.read_table(config["subset"], index_col=["sample"], dtype=str).sort_index(level=0)
    idr_list = []
    idr_list = idr_list + expand("idr/distal/{unit.sample}/{unit.sample}-peakset-distal.narrowPeak", unit=units[units['unit'].notnull()].reset_index().drop_duplicates(subset='sample').itertuples())
    #idr_list = idr_list + expand("idr/distal/{unit.sample}/{unit.sample}-peakset-distal.narrowPeak", unit=units[units['unit'].isnull()].reset_index().drop_duplicates(subset='sample').itertuples())
    return(idr_list)

def get_idr_subset_input_promoter(wildcards):
    units = pd.read_table(config["subset"], index_col=["sample"], dtype=str).sort_index(level=0)
    idr_list = []
    idr_list = idr_list + expand("idr/promoter/{unit.sample}/{unit.sample}-peakset-promoter.narrowPeak", unit=units[units['unit'].notnull()].reset_index().drop_duplicates(subset='sample').itertuples())
    #idr_list = idr_list + expand("idr/promoter/{unit.sample}/{unit.sample}-peakset-promoter.narrowPeak", unit=units[units['unit'].isnull()].reset_index().drop_duplicates(subset='sample').itertuples())
    return(idr_list)


rule merge_macs_peaks:
    input:
        get_idr_input
    output:
        "macs/bedgraph/{sample}_peaks.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} | bedtools merge -i - -c 5 -o mean > {output}
        """

rule unionbedgraph:
    input:
        expand("macs/bedgraph/{unit.sample}_peaks.bedgraph", unit=subset.reset_index().drop_duplicates(subset='sample').itertuples())
    output:
        "macs/unionbedg/union.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    params:
        units=units
    shell:
        """
        bedtools unionbedg -i {input} > {output}.unsorted && sort -k1,1 -k2,2n {output}.unsorted | bedtools merge -i - > {output}
        """


rule unionbedgraph_TSS:
    input:
        expand("macs/bedgraph/{unit.sample}_peaks.bedgraph", unit=units.reset_index().drop_duplicates(subset='sample').itertuples())
    output:
        "macs/unionbedg/union_TSS.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    params:
        units=units,
        TSS="/mnt/data1/John/DATA_Common/mm10_refseq.bed"
    shell:
        """
        bedtools unionbedg -i {input} > {output}.unsorted && sort -k1,1 -k2,2n {output}.unsorted | bedtools merge -i - | bedtools intersect -a - -b {params.TSS} -wa -u > {output}
        """

rule coverage_unit:
    input:
        bed="macs/unionbedg/union.bedgraph",
        bam="bam/units/{sample}-{unit}.bam"
    output:
        "macs/coverage/{sample}-{unit}-peaks_coverage.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    params: ""
    log: "logs/bedtools/coverage/{sample}-{unit}.log"
    shell:
        """
        samtools index {input.bam} && bedtools multicov -bams {input.bam} -bed {input.bed} 2> {log} 1> {output}
        """

rule coverage_sample:
    input:
        bed="macs/unionbedg/union.bedgraph",
        bam=get_bam_input
    output:
        "macs/coverage/{sample}-peaks_coverage.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    params: ""
    log: "logs/bedtools/coverage/{sample}.log"
    shell:
        """
        samtools index {input.bam} && bedtools multicov -bams {input.bam} -bed {input.bed} 2> {log} 1> {output}
        """

rule coverage_unit_TSS:
    input:
        bed="macs/unionbedg/union_TSS.bedgraph",
        bam="bam/units/{sample}-{unit}.bam"
    output:
        "macs/coverage/{sample}-{unit}-TSS_coverage.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    params: ""
    shell:
        """
        samtools index {input.bam} && bedtools multicov -bams {input.bam} -bed {input.bed} > {output}
        """

rule coverage_sample_TSS:
    input:
        bed="macs/unionbedg/union_TSS.bedgraph",
        bam=get_bam_input
    output:
        "macs/coverage/{sample}-TSS_coverage.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    params: ""
    shell:
        """
        samtools index {input.bam} && bedtools multicov -bams {input.bam} -bed {input.bed} > {output}
        """

