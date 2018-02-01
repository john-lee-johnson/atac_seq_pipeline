rule deseq2_to_bed:
    input:
        "results/diffexp/{contrast}_{result}.diffexp.tsv",
    output:
        "results/diffexp/bed/{contrast}_{result}.diffexp.bed",
    shell:
        """
        tail -n +2 {input} | cut -d' ' -f1 | cut -d'.' -f1,2,3 --output-delimiter=$'\t' > {output}
        """

rule homer:
    input:
        "results/diffexp/bed/{contrast}_{result}.diffexp.bed",
    output:
        "motif/{contrast}_{result}_sites/homerResults.html",
    conda:
        "../envs/homer.yaml"
    params:
        dir="motif/{contrast}_{result}_sites/",
        genome=config["params"]["homer"]
    threads: 24
    shell:
        """
        findMotifsGenome.pl {input} {params.genome} {params.dir} -p {threads}
        """
