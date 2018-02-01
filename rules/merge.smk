def get_unit_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    unit_list = units.loc[[wildcards.sample]]['unit'].values
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        unit_list = expand("bam/units/{sample}-{unit}.bam", unit=unit_list, sample=wildcards.sample)
    return(unit_list)

def get_lane_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[wildcards.sample]['unit'].notnull().any():
        units = units.reset_index().set_index(['sample', 'unit'])
        if units.loc[(wildcards.sample, wildcards.unit)][['lane']].shape[0] >1:
            lane_list = units.loc[(wildcards.sample, wildcards.unit)]['lane'].values
            for lane in lane_list:
                bam_list.append("star/{sample}-{unit}-{lane}/Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample, unit=wildcards.unit, lane = lane))
    else:
        lane_list = units.loc[wildcards.sample]['lane'].values
        for lane in lane_list:
            bam_list.append("star/{sample}-{lane}/Aligned.sortedByCoord.out.bam".format(sample=wildcards.sample, lane = lane))
    return bam_list

rule merge_unit_lane:
    input:
        get_lane_input
    output:
        "merge/{sample}-{unit}.bam"
    log:
        "logs/picard/mergebam/{sample}-{unit}.log"
    params: ""
    wrapper:
        "master/bio/picard/mergesamfiles"

rule merge_sample_lane:
    input:
        get_lane_input
    output:
        "merge/{sample}.bam"
    log:
        "logs/picard/mergebam/{sample}.log"
    params: ""
    wrapper:
        "master/bio/picard/mergesamfiles"


rule merge_unit:
    input:
        get_unit_input
    output:
        "bam/merge/{sample}.bam"
    log:
        "logs/picard/mergebam/{sample}.log"
    params: ""
    wrapper:
        "master/bio/picard/mergesamfiles"

rule merge_all:
    input:
        expand("bam/{unit.sample}.bam", unit=units.reset_index().itertuples())
    output:
        "bam/all_samples_merged.bam"
    log:
        "logs/picard/mergebam/all_samples_merged.log"
    params: ""
    wrapper:
        "master/bio/picard/mergesamfiles"


