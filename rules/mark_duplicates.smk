def get_units(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if (units.loc[[wildcards.sample]]['unit'].notnull().any() and hasattr(wildcards, 'unit')):
    # has unit
        units = units.reset_index().set_index(['sample', 'unit'])
        if units.loc[(wildcards.sample, wildcards.unit)][['lane']].shape[0] >1:
            # unit has multiple fastq
            bam_list.append("merge/{sample}-{unit}.bam".format(sample=wildcards.sample, unit=wildcards.unit))
        else:
            # unit has one fastq
            bam_list.append("star/{sample}-{unit}-{lane}/Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample, unit=wildcards.unit, lane=units.loc[[(wildcards.sample, wildcards.unit)]]['lane'].values[0]))
    elif units.loc[[wildcards.sample]]['unit'].isnull().any():
        # no unit
        if units.loc[[wildcards.sample]][['lane']].shape[0] >1:
            # sample has multiple fastq
            bam_list.append("merge/{sample}.bam".format(sample=wildcards.sample))
        else:
            # sample has one fastq
            bam_list.append("star/{sample}-{lane}/Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample, lane=units.loc[[wildcards.sample]]['lane'].values[0]))
    #print(bam_list)
    return bam_list

rule mark_duplicates_units:
    input:
        get_units
    output:
        bam=temp("dedup/unfiltered/{sample}-{unit}.unfiltered.bam"),
        metrics="dedup/unfiltered/{sample}-{unit}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}-{unit}.log"
    params:
        "REMOVE_DUPLICATES=true "
        "READ_NAME_REGEX=null "
        "TMP_DIR=/mnt/data1/John/tmp "
    wrapper:
        "master/bio/picard/markduplicates"

rule mark_duplicates_samples:
    input:
        get_units
    output:
        bam=temp("dedup/unfiltered/{sample}.unfiltered.bam"),
        metrics="dedup/unfiltered/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        "REMOVE_DUPLICATES=true "
        "READ_NAME_REGEX=null "
        "TMP_DIR=/mnt/data1/John/tmp "
    wrapper:
        "master/bio/picard/markduplicates"

