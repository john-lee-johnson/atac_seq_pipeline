def get_sample_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        bam_list.append("bam/merge/{sample}.bam".format(sample=wildcards.sample))
    else:
        bam_list.append("bam/{sample}.bam".format(sample=wildcards.sample))
    return bam_list


rule insert_size_sample:
    input:
        get_sample_input
    output:
        txt="stats/{sample}.isize.txt",
        pdf="stats/{sample}.isize.pdf"
    log:
        "logs/picard/insert_size/{sample}.log"
    conda:
        "../../../github/snakemake-wrappers/bio/picard/collectinsertsizemetrics/environment.yaml"
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE "
        "TMP_DIR=tmp "
    wrapper:
        "file://../../github/snakemake-wrappers/bio/picard/collectinsertsizemetrics/"

rule insert_size_unit:
    input:
        "bam/units/{sample}-{unit}.bam"
    output:
        txt="stats/{sample}-{unit}.isize.txt",
        pdf="stats/{sample}-{unit}.isize.pdf"
    log:
        "logs/picard/insert_size/{sample}-{unit}.log"
    conda:
        "../../../github/snakemake-wrappers/bio/picard/collectinsertsizemetrics/environment.yaml"
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE "
        "TMP_DIR=tmp "
    wrapper:
        "file://../../github/snakemake-wrappers/bio/picard/collectinsertsizemetrics/"

rule insert_size_all:
    input:
        "merge/All_Samples_Merged.bam"
    output:
        txt="stats/All_Samples_Merged.isize.txt",
        pdf="stats/All_Samples_Merged.isize.pdf"
    log:
        "logs/picard/insert_size/All_Samples_Merged.log"
    conda:
        "../../../github/snakemake-wrappers/bio/picard/collectinsertsizemetrics/environment.yaml"
    params:
        # optional parameters (e.g. relax checks as below)
        "VALIDATION_STRINGENCY=LENIENT "
        "METRIC_ACCUMULATION_LEVEL=null "
        "METRIC_ACCUMULATION_LEVEL=SAMPLE "
        "TMP_DIR=tmp "
    wrapper:
        "file://../../github/snakemake-wrappers/bio/picard/collectinsertsizemetrics/"
