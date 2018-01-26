def get_rep_info(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    unit_list = units.loc[[wildcards.sample]]['unit'].unique()
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        units_list = expand("macs/units/{sample}-{unit}_peaks.narrowPeak", unit=unit_list, sample=wildcards.sample)
        units_list.append("macs/{wildcards.sample}_peaks.narrowPeak".format(wildcards=wildcards))
    return(units_list)

rule idr_optimal:
    input:
        "macs/{sample}_peaks.narrowPeak", "macs/pseudoreps/{sample}-pseurep1_peaks.narrowPeak", "macs/pseudoreps/{sample}-pseurep2_peaks.narrowPeak"
    output:
        "idr/{sample}/{sample}-optimal-peakset.narrowPeak"
    log:
        "logs/idr/{sample}-optimal.log"
    conda:
        "../envs/idr.yaml"
    params: ""
    shell:
        """
        idr -s {input[1]} {input[2]} --peak-list {input[0]} --input-file-type narrowPeak --rank signal.value --output-file {output} --output-file-type narrowPeak --log-output-file {log} --plot --use-best-multisummit-IDR --idr-threshold 0.05
        """

rule idr_conservative:
    input:
        get_rep_info
    output:
        "idr/{sample}/{sample}-conservative-peakset.narrowPeak"
    log:
        "logs/idr/{sample}-conservative.log"
    conda:
        "../envs/idr.yaml"
    params: ""
    shell:
        """
        idr -s {input[0]} {input[1]} --peak-list {input[2]} --input-file-type narrowPeak --rank signal.value --output-file {output} --output-file-type narrowPeak --log-output-file {log} --plot --use-best-multisummit-IDR --idr-threshold 0.05
        """

