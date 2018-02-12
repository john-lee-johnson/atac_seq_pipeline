def get_sample_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        # has units
        bam_list.append("bam/merge/{sample}.bam".format(sample=wildcards.sample))
    else:
        # no units
        bam_list.append("bam/{sample}.bam".format(sample=wildcards.sample))
    return bam_list

def get_bdg_input(wildcards):
    units = pd.read_table(config["units"], index_col=["sample"], dtype=str).sort_index(level=0)
    bam_list = []
    if units.loc[[wildcards.sample]]['unit'].notnull().values.any():
        # has units
        bam_list.append("units/{sample}-{unit}_treat_pileup.bdg".format(sample=wildcards.sample))
    else:
        # no units
        bam_list.append("dedup/{sample}.bam".format(sample=wildcards.sample))
    return bam_list

rule call_peaks_sample:
    input: get_sample_input
    output:
        peak = "macs/{sample}_peaks.narrowPeak",
        bdg = "macs/{sample}_treat_pileup.bdg",
        ctrl = "macs/{sample}_control_lambda.bdg"
    log:
        "logs/macs2/call_peaks/{sample}.log"
    conda:
        "../envs/macs2.yaml"
    params:
        name = "{sample}",
        extra = "{}".format(config["params"]["macs2"])
    shell:
        """
        macs2 callpeak -t {input} -f BAMPE -g {params.extra} --outdir macs -n {params.name} --nomodel --keep-dup all -B --call-summits -q 0.9 &> {log}
        """

rule call_peaks_units:
    input: "bam/units/{sample}-{unit}.bam"
    output:
        peak = "macs/units/{sample}-{unit}_peaks.narrowPeak",
        bdg = "macs/units/{sample}-{unit}_treat_pileup.bdg",
        ctrl = "macs/units/{sample}-{unit}_control_lambda.bdg"
    log: "logs/macs2/call_peaks/{sample}-{unit}.log"
    conda: "../envs/macs2.yaml"
    params:
        name = "{sample}-{unit}",
        extra = "{}".format(config["params"]["macs2"])
    shell:
        """
        macs2 callpeak -t {input} -f BAMPE -g {params.extra} --outdir macs/units -n {params.name} --nomodel --keep-dup all -B --call-summits -q 0.9 &> {log}
        """

rule call_peaks_pseudoreps:
    input: "dedup/split/{sample}-{pseurep}.bed"
    output:
        peak = temp("macs/pseudoreps/{sample}-{pseurep}_peaks.narrowPeak"),
        bdg = temp("macs/pseudoreps/{sample}-{pseurep}_treat_pileup.bdg"),
        ctrl = temp("macs/pseudoreps/{sample}-{pseurep}_control_lambda.bdg"),
        summits = temp("macs/pseudoreps/{sample}-{pseurep}_summits.bed"),
        xls = temp("macs/pseudoreps/{sample}-{pseurep}_peaks.xls")
    log: "logs/macs2/call_peaks/{sample}-{pseurep}.log"
    conda: "../envs/macs2.yaml"
    params:
        name = "{sample}-{pseurep}",
        extra = "{}".format(config["params"]["macs2"])
    shell:
        """
        macs2 callpeak -t {input} -f BEDPE -g {params.extra} --outdir macs/pseudoreps -n {params.name} --nomodel --keep-dup all -B --call-summits -q 0.9 &> {log}
        """

rule macs_bigwig_sample:
    input: "macs/{sample}_treat_pileup.bdg"
    output: "bw/macs2/{sample}_treat_pileup.bw"
    log:
        "logs/macs2/bw/{sample}.log"
    conda:
        "../envs/macs2.yaml"
    params:
        genome = "{}".format(config["params"]["genome"])
    log:
        "logs/macs2/bw/{sample}.log"
    shell:
        """
        LC_COLLATE=C sort -k 1,1 -k 2,2n {input} > {input}.sorted && bedGraphToBigWig {input}.sorted {params.genome} {output} && rm {input}.sorted &> {log}
        """

rule macs_bigwig_unit:
    input: "macs/units/{sample}-{unit}_treat_pileup.bdg"
    output: "bw/macs2/units/{sample}-{unit}_treat_pileup.bw"
    log: "logs/macs2/bw/{sample}-{unit}.log"
    conda: "../envs/macs2.yaml"
    params:
        genome = "{}".format(config["params"]["genome"])
    log:
        "logs/macs2/bw/{sample}-{unit}.log"
    shell:
        """
        LC_COLLATE=C sort -k1,1 -k2,2n {input} > {input}.sorted && bedGraphToBigWig {input}.sorted {params.genome} {output} && rm {input}.sorted &> {log}
        """

