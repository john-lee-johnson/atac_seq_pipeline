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

rule merge_idr_peaks_distal:
    input:
        "idr/distal/{sample}/{sample}-peakset-distal.narrowPeak"
    output:
        "idr/bedgraph/distal/{sample}_peaks.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} | bedtools merge -i - -c 5 -o mean > {output}
        """

rule merge_idr_peaks_promoter:
    input:
        "idr/promoter/{sample}/{sample}-peakset-promoter.narrowPeak"
    output:
        "idr/bedgraph/promoter/{sample}_peaks.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        sort -k1,1 -k2,2n {input} | bedtools merge -i - -c 5 -o mean > {output}
        """

rule unionbedgraph_distal:
    input:
        expand("idr/bedgraph/distal/{unit.sample}_peaks.bedgraph", unit=subset.reset_index().drop_duplicates(subset='sample').itertuples())
    output:
        "idr/union/union-distal.bedgraph"
    conda:
        "../envs/bedtools.yaml"
    params:
        units=units
    shell:
        """
        bedtools unionbedg -i {input} > {output}.unsorted && sort -k1,1 -k2,2n {output}.unsorted | bedtools merge -i - > {output}
        """


rule unionbedgraph_promoter:
    input:
        expand("idr/bedgraph/promoter/{unit.sample}_peaks.bedgraph", unit=subset.reset_index().drop_duplicates(subset='sample').itertuples())
    output:
        "idr/union/union-promoter.bedgraph",
    conda:
        "../envs/bedtools.yaml"
    params:
        units=units
    shell:
        """
        bedtools unionbedg -i {input} > {output}.unsorted && sort -k1,1 -k2,2n {output}.unsorted | bedtools merge -i - > {output}
        """


rule idr_qvalue_union_distal:
    input:
        union="idr/union/union-distal.bedgraph",
        idr=expand("idr/distal/{unit.sample}/{unit.sample}-peakset-distal.narrowPeak", unit=subset[subset['unit'].notnull()].reset_index().drop_duplicates(subset='sample').itertuples())
    output:
        "idr/union/union-fdr-distal.bed"
    conda:
        "../envs/bedtools.yaml"
    params: ""
    log: "logs/union/union-fdr.log"
    script:
        "../scripts/union-fdr.py"


rule idr_qvalue_union_promoter:
    input:
        union="idr/union/union-promoter.bedgraph",
        idr=expand("idr/distal/{unit.sample}/{unit.sample}-peakset-promoter.narrowPeak", unit=subset[subset['unit'].notnull()].reset_index().drop_duplicates(subset='sample').itertuples())
    output:
        "idr/union/union-fdr-promoter.bed"
    conda:
        "../envs/bedtools.yaml"
    params: ""
    log: "logs/union/union-fdr.log"
    script:
        "../scripts/union-fdr.py"


rule heatmap_distal:
    input:
        "idr/union/union-fdr-distal.bed"
    output:
        output1 = "plots/heatmap_distal.png",
        outputk1 = "plots/heatmap_distal_kcenters.png",
        output2 = "plots/heatmap2_distal.png",
        outputk2 = "plots/heatmap2_distal_kcenters.png"
    conda:
        "../envs/deseq2.yaml"
    params:
        k1 = config["heatmap-distal"]["k1"],
        k2 = config["heatmap-distal"]["k2"],
        label = config["heatmap"]["label"],
        order = config["heatmap"]["order"],
        remove = config["heatmap-distal"]["remove"],
        order2 = config["heatmap-distal"]["order2"],
    log: "logs/union/heatmap-distal.log"
    script:
        "../scripts/heatmap.R"


rule heatmap_proximal:
    input:
        "idr/union/union-fdr-promoter.bed"
    output:
        output1 = "plots/heatmap_promoter.png",
        outputk1 = "plots/heatmap_promoter_kcenters.png",
        output2 = "plots/heatmap2_promoter.png",
        outputk2 = "plots/heatmap2_promoter_kcenters.png"
    conda:
        "../envs/deseq2.yaml"
    params:
        k1 = config["heatmap-promoter"]["k1"],
        k2 = config["heatmap-promoter"]["k2"],
        label = config["heatmap"]["label"],
        order = config["heatmap"]["order"],
        remove = config["heatmap-promoter"]["remove"],
        order2 = config["heatmap-promoter"]["order2"],
    log: "logs/union/heatmap-promoter.log"
    script:
        "../scripts/heatmap.R"

