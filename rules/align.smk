def get_trimmed(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    if not is_single_end(wildcards.sample, wildcards.lane):
        # paired-end sample
        if (units.loc[[wildcards.sample]]['unit'].notnull().values.any() and hasattr(wildcards, 'unit')):
           # has units
            return expand("trimmed/{sample}-{unit}-{lane}_{group}.fastq.gz",
                      group=[1,2], **wildcards)
        else:
            # no units
            return expand("trimmed/{sample}-{lane}_{group}.fastq.gz",
                      group=[1,2], **wildcards)
    # single end sample
    else:
        if (units.loc[[wildcards.sample]]['unit'].notnull().values.any() and hasattr(wildcards, 'unit')):
            # has units
            return "trimmed/{sample}-{unit}-{lane}.fastq.gz".format(**wildcards)
            # no units
    return "trimmed/{sample}-{lane}.fastq.gz".format(**wildcards)

rule align_no_unit:
    input:
        sample=get_trimmed
    output:
        # see STAR manual for additional output files
        protected("star/{sample}-{lane}/Aligned.sortedByCoord.out.bam"),
    log:
        "logs/star/{sample}-{lane}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="{}".format(config["params"]["star"])
    threads: 17
    wrapper:
        "https://bitbucket.org/john-lee-johnson/snakemake-wrappers/raw/master/bio/star/align"

rule align_unit:
    input:
        sample=get_trimmed
    output:
        # see STAR manual for additional output files
        protected("star/{sample}-{unit}-{lane}/Aligned.sortedByCoord.out.bam"),
    log:
        "logs/star/{sample}-{unit}-{lane}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="{}".format(config["params"]["star"])
    threads: 17
    wrapper:
        "https://bitbucket.org/john-lee-johnson/snakemake-wrappers/raw/master/bio/star/align"
