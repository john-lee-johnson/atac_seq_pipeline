def get_sample_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        # has units
        bam_list.append("bam/merge/{sample}.bed".format(sample=wildcards.sample))
    else:
        # no units
        bam_list.append("bam/{sample}.bed".format(sample=wildcards.sample))
    return bam_list


rule pseudoreplicates:
    input:
        "bam/bed/{sample}.bed"
    output:
        temp("dedup/split/{sample}-pseurep1.bed"),
        temp("dedup/split/{sample}-pseurep2.bed")
    params:
        name = "bed/{sample}",
        outname = "dedup/split/{sample}-pseurep"
    shell:
        """
        shuf {input} > {input}.shuffle && cut -d $'\t' -f 1,2,6 {input}.shuffle > {input}.shuffle.bed && split --number=l/2 {input}.shuffle.bed {params.outname} --additional-suffix=.bed --numeric-suffixes=1 --suffix-length=1 && rm {input}.shuffle.bed && rm {input}.shuffle
        """
