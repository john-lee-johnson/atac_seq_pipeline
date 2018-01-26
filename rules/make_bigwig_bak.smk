rule make_bigwig_unit:
    input:
        "dedup/{sample}-{unit}.bam"
    output:
        "bw/{sample}-{unit}.bw"
    log:
        "logs/deeptools/{sample}-{unit}.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--normalizeUsingRPKM",
        blacklist = "/mnt/data1/John/DATA_Common/blacklist/mm10.blacklist.bed"
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} --blackListFileName {params.blacklist} 2> {log}
        """

rule make_bigwig_sample:
    input:
        "dedup/{sample}.bam"
    output:
        "bw/{sample}.bw"
    log:
        "logs/deeptools/{sample}.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--normalizeUsingRPKM",
        blacklist = "/mnt/data1/John/DATA_Common/blacklist/mm10.blacklist.bed"
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} --blackListFileName {params.blacklist} 2> {log}
        """
