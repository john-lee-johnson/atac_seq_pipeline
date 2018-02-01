def get_scaling_factor(wildcards):
    factors = pd.read_table("counts/sizeFactors_allPeaks.tsv", index_col=0, sep=' ', header = None, names = ["sample", "scaleFactor"])
    name = wildcards.sample
    factors['scaleFactor'] = factors['scaleFactor'].pow(-1)
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        name += "-"
        name += wildcards.unit
        return factors.loc[name]['scaleFactor']
    elif units.loc[[wildcards.sample]]['unit'].isnull().values.any():
        return factors.loc[name]['scaleFactor']

def get_scaling_factor_TSS(wildcards):
    factors = pd.read_table("counts/sizeFactors_TSS.tsv", index_col=0, sep=' ', header = None, names = ["sample", "scaleFactor"])
    name = wildcards.sample
    factors['scaleFactor'] = factors['scaleFactor'].pow(-1)
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        name += "-"
        name += wildcards.unit
        return factors.loc[name]['scaleFactor']
    elif units.loc[[wildcards.sample]]['unit'].isnull().values.any():
        return factors.loc[name]['scaleFactor']


def get_bam_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        bam_list.append("bam/merge/{sample}.bam".format(sample=wildcards.sample))
    else:
        bam_list.append("bam/{sample}.bam".format(sample=wildcards.sample))
    return bam_list


rule make_bigwig_unit:
    input:
        "bam/units/{sample}-{unit}.bam"
    output:
        "bw/{sample}-{unit}-scaleFactor.bw"
    log:
        "logs/deeptools/{sample}-{unit}-scaleFactor.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--scaleFactor",
        scale = get_scaling_factor,
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} {params.scale} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """

rule make_bigwig_sample:
    input:
        get_bam_input
    output:
        "bw/{sample}-scaleFactor.bw"
    log:
        "logs/deeptools/{sample}-scaleFactor.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--scaleFactor",
        scale = get_scaling_factor,
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} {params.scale} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """

rule make_bigwig_unit_TSS:
    input:
        bam = "bam/units/{sample}-{unit}.bam",
        scale = "counts/sizeFactors_TSS.tsv"
    output:
        "bw/{sample}-{unit}-scaleFactor_TSS.bw"
    log:
        "logs/deeptools/{sample}-{unit}-scaleFactor_TSS.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--scaleFactor",
        scale = get_scaling_factor_TSS,
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input.bam} && bamCoverage --bam {input.bam} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} {params.scale} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """

rule make_bigwig_sample_TSS:
    input:
        bam = get_bam_input,
        scale = "counts/sizeFactors_TSS.tsv"
    output:
        "bw/{sample}-scaleFactor_TSS.bw"
    log:
        "logs/deeptools/{sample}-scaleFactor_TSS.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--scaleFactor",
        scale = get_scaling_factor_TSS,
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input.bam} && bamCoverage --bam {input.bam} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} {params.scale} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """


rule make_bigwig_unit_rpkm:
    input:
        "bam/units/{sample}-{unit}.bam"
    output:
        "bw/rpkm/units/{sample}-{unit}-RPKM.bw"
    log:
        "logs/deeptools/{sample}-{unit}-RPKM.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--normalizeUsingRPKM",
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """

rule make_bigwig_sample_rpkm:
    input:
        get_bam_input
    output:
        "bw/rpkm/{sample}-RPKM.bw"
    log:
        "logs/deeptools/{sample}-RPKM.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--normalizeUsingRPKM",
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """


rule make_bedgraph_unit:
    input:
        "bam/units/{sample}-{unit}.bam"
    output:
        "bw/{sample}-{unit}.bedGraph"
    log:
        "logs/deeptools/{sample}-{unit}.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bedgraph",
        normalize = "--scaleFactor",
        scale = get_scaling_factor,
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} {params.scale} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """

rule make_bedgraph_sample:
    input:
        get_bam_input
    output:
        "bw/{sample}.bedGraph"
    log:
        "logs/deeptools/{sample}.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bedgraph",
        normalize = "--scaleFactor",
        scale = get_scaling_factor,
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} {params.scale} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """


rule make_bigwig_all:
    input:
        "bam/all_samples_merged.bam"
    output:
        "bw/all_samples_merged.bw"
    log:
        "logs/deeptools/all_samples_merged.log"
    conda:
        "../envs/deeptools.yaml"
    params:
        binsize = 10,
        outformat = "bigwig",
        normalize = "--normalizeUsingRPKM",
        blacklist = config["params"]["blacklist"],
    threads: 1
    shell:
        """
        samtools index {input} && bamCoverage --bam {input} -o {output} --binSize {params.binsize} --outFileFormat {params.outformat} {params.normalize} --blackListFileName {params.blacklist} -p {threads} 2> {log}
        """
