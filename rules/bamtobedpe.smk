def get_bam_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
    # has units
        return "bam/merge/{sample}.bam".format(sample=wildcards.sample)
    else:
    # no units
        return "bam/{sample}.bam".format(sample=wildcards.sample)

rule bamtobedpe:
    input:
        get_bam_input
    output:
        temp("bam/bed/{sample}.bed")
    params: ""
    conda:
        "../envs/bedtools.yaml"
    shell:
        """
        samtools sort -n {input} -T /mnt/data1/John/tmp | samtools view -bf 0x2 > {input}.sorted && samtools fixmate {input}.sorted {input}.fixed && bedtools bamtobed -i {input}.fixed -bedpe > {output} && rm {input}.sorted && rm {input}.fixed
        """

