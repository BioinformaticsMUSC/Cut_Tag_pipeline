"""
Quality Control Rules
"""

rule fastqc:
    """Run FastQC on raw reads"""
    input:
        "data/{sample}_R{read}.fastq.gz"
    output:
        html="results/fastqc/{sample}_R{read}_fastqc.html",
        zip="results/fastqc/{sample}_R{read}_fastqc.zip"
    params:
        outdir="results/fastqc"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        "fastqc -t {threads} -o {params.outdir} {input}"

rule fastqc_trimmed:
    """Run FastQC on trimmed reads"""
    input:
        "results/trimmed/{sample}_R{read}_trimmed.fastq.gz"
    output:
        html="results/fastqc_trimmed/{sample}_R{read}_trimmed_fastqc.html",
        zip="results/fastqc_trimmed/{sample}_R{read}_trimmed_fastqc.zip"
    params:
        outdir="results/fastqc_trimmed"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        "fastqc -t {threads} -o {params.outdir} {input}"

rule multiqc:
    """Generate MultiQC report"""
    input:
        expand("results/fastqc/{sample}_R{read}_fastqc.zip", sample=SAMPLES, read=[1,2]),
        expand("results/fastqc_trimmed/{sample}_R{read}_trimmed_fastqc.zip", sample=SAMPLES, read=[1,2]),
        expand("results/aligned/{sample}.sorted.bam", sample=SAMPLES)
    output:
        "results/qc/multiqc_report.html"
    params:
        outdir="results/qc",
        title="CUT&Tag Analysis Report"
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc -o {params.outdir} -n multiqc_report.html "
        "-i '{params.title}' results/fastqc results/fastqc_trimmed results/aligned"