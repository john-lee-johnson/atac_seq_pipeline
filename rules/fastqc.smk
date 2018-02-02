rule fastqc:
    input:
        "fastq/{fastq}.fastq.gz"
    output:
        html="qc/{fastq}_fastqc.html",
        zip="qc/{fastq}_fastqc.zip",
    params: ""
    wrapper:
        "master/bio/fastqc"

rule fastqc_trimmed_with_unit:
    input:
        "trimmed/{sample}-{unit}-{fastq}.fastq.gz"
    output:
        html="qc/{sample}-{unit}-{fastq}_trimmed_fastqc.html",
        zip="qc/{sample}-{unit}-{fastq}_trimmed_fastqc.zip"
    params: ""
    wrapper:
        "master/bio/fastqc"

rule fastqc_trimmed_without_unit:
    input:
        "trimmed/{sample}-{fastq}.fastq.gz"
    output:
        html="qc/{sample}-{fastq}_trimmed_fastqc.html",
        zip="qc/{sample}-{fastq}_trimmed_fastqc.zip"
    params: ""
    wrapper:
        "master/bio/fastqc"
