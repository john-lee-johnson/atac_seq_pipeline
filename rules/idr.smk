def get_rep_info(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    unit_list = units.loc[[wildcards.sample]]['unit'].unique()
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        units_list = expand("macs/units/{sample}-{unit}_peaks.narrowPeak", unit=unit_list, sample=wildcards.sample)
        units_list.append("macs/{wildcards.sample}_peaks.narrowPeak".format(wildcards=wildcards))
    return(units_list)

def get_idr_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    #unit_list = units.loc[[wildcards.sample]]['unit'].unique()
    units_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        units_list.append("idr/{sample}/{sample}-conservative-peakset.narrowPeak".format(sample=wildcards.sample))
    else:
        units_list.append("idr/{sample}/{sample}-optimal-peakset.narrowPeak".format(sample=wildcards.sample))
    return(units_list)


rule idr_optimal:
    input:
        "macs/{sample}_peaks.narrowPeak", "macs/pseudoreps/{sample}-pseurep1_peaks.narrowPeak", "macs/pseudoreps/{sample}-pseurep2_peaks.narrowPeak"
    output:
        peak="idr/{sample}/{sample}-optimal-peakset.narrowPeak",
        tabix="idr/{sample}/{sample}-optimal-peakset.narrowPeak.gz",
    log:
        "logs/idr/{sample}-optimal.log"
    conda:
        "../envs/idr.yaml"
    params: ""
    shell:
        """
        idr -s {input[1]} {input[2]} --peak-list {input[0]} --input-file-type narrowPeak --rank signal.value --output-file {output.peak} --output-file-type narrowPeak --log-output-file {log} --plot --use-best-multisummit-IDR --idr-threshold 0.05 && sort {output.peak} -k1,1 -k2,2n | bgzip > {output.tabix} && tabix -p bed {output.tabix}

        """

rule idr_conservative:
    input:
        get_rep_info
    output:
        peak="idr/{sample}/{sample}-conservative-peakset.narrowPeak",
        tabix="idr/{sample}/{sample}-conservative-peakset.narrowPeak.gz"
    log:
        "logs/idr/{sample}-conservative.log"
    conda:
        "../envs/idr.yaml"
    params: ""
    shell:
        """
        idr -s {input[0]} {input[1]} --peak-list {input[2]} --input-file-type narrowPeak --rank signal.value --output-file {output.peak} --output-file-type narrowPeak --log-output-file {log} --plot --use-best-multisummit-IDR --idr-threshold 0.05 && sort {output.peak} -k1,1 -k2,2n | bgzip > {output.tabix} && tabix -p bed {output.tabix}
        """

rule split_distal_peaks:
    input:
        get_idr_input
    output:
        distal="idr/distal/{sample}/{sample}-peakset-distal.narrowPeak",
        promoter="idr/promoter/{sample}/{sample}-peakset-promoter.narrowPeak"
    log:
        "logs/idr/{sample}-split-distal.log"
    conda:
        "../envs/bedtools.yaml"
    params: "/mnt/data1/John/DATA_Common/GENCODE/M15_GRCm38.p5/Gencode_Gene_Collapsed_to_Single_Entry_Promoter_Only_Sorted.bed"
    shell:
        """
        bedtools intersect -a {input} -b {params} -v -nonamecheck > {output.distal} && bedtools intersect -a {input} -b {params} -wa -nonamecheck > {output.promoter}
        """
