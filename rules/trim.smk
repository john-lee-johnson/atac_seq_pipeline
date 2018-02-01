def get_fastq(wildcards):
    units = pd.read_table(config["units"], index_col=["sample", "lane"], dtype=str).sort_index(level=0)
    return "fastq/"+units.loc[(wildcards.sample, wildcards.lane), ["fq1", "fq2"]].dropna()

rule cutadapt_pe_with_unit:
    input:
        get_fastq
    output:
        fastq1=temp("trimmed/{sample}-{unit}-{lane}_1.fastq.gz"),
        fastq2=temp("trimmed/{sample}-{unit}-{lane}_2.fastq.gz"),
        qc="trimmed/{sample}-{unit}-{lane}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-pe"])
    log:
        "logs/cutadapt/{sample}-{unit}-{lane}.log"
    wrapper:
        "master/bio/cutadapt/pe"

rule cutadapt_pe_no_unit:
    input:
        get_fastq
    output:
        fastq1=temp("trimmed/{sample}-{lane}_1.fastq.gz"),
        fastq2=temp("trimmed/{sample}-{lane}_2.fastq.gz"),
        qc="trimmed/{sample}-{lane}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-pe"])
    log:
        "logs/cutadapt/{sample}-{lane}.log"
    wrapper:
        "master/bio/cutadapt/pe"


rule cutadapt_with_unit:
    input:
        get_fastq
    output:
        fastq=temp("trimmed/{sample}-{unit}-{lane}.fastq.gz"),
        qc="trimmed/{sample}-{unit}-{lane}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-se"])
    log:
        "logs/cutadapt/{sample}-{unit}-{lane}.log"
    wrapper:
        "master/bio/cutadapt/se"

rule cutadapt_no_unit:
    input:
        get_fastq
    output:
        fastq=temp("trimmed/{sample}-{lane}.fastq.gz"),
        qc="trimmed/{sample}-{lane}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-se"])
    log:
        "logs/cutadapt/{sample}-{lane}.log"
    wrapper:
        "master/bio/cutadapt/se"
