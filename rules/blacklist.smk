def get_units(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].isnull().values.any():
        bam_list.append("dedup/unfiltered/{sample}.unfiltered.bam".format(sample=wildcards.sample))
    return bam_list


rule remove_blacklist_unit:
    input:
        "dedup/unfiltered/{sample}-{unit}.unfiltered.bam"
    output:
        cleaned="bam/units/{sample}-{unit}.bam",
        failed="dedup/failed/{sample}-{unit}.blacklist.failed.reads"
    log:
        "logs/bamutils/{sample}-{unit}.log"
    conda:
        "../envs/bamutils.yaml"
    params:
        "{}".format(config["params"]["blacklist"])
    shell:
        """
        samtools index {input} && bamutils filter {input} {output.cleaned} -excludebed {params} nostrand -failed {output.failed}
        """

rule remove_blacklist_sample:
    input:
        get_units
    output:
        cleaned="bam/{sample}.bam",
        failed="dedup/failed/{sample}.blacklist.failed.reads"
    log:
        "logs/bamutils/{sample}.log"
    conda:
        "../envs/bamutils.yaml"
    params:
        "{}".format(config["params"]["blacklist"])
    shell:
        """
        samtools index {input} && bamutils filter {input} {output.cleaned} -excludebed {params} nostrand -failed {output.failed}
        """
