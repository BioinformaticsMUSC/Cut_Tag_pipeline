"""
Read Trimming Rules
"""

rule trimmomatic:
    """Trim adapters and low-quality bases from paired-end reads"""
    input:
        r1=lambda wildcards: f"data/{wildcards.sample}_R1{config['input']['fastq_suffix']}",
        r2=lambda wildcards: f"data/{wildcards.sample}_R2{config['input']['fastq_suffix']}",
        adapters=config["trimmomatic"]["adapters"]
    output:
        r1="results/trimmed/{sample}_R1_trimmed.fastq.gz",
        r2="results/trimmed/{sample}_R2_trimmed.fastq.gz",
        r1_unpaired="results/trimmed/{sample}_R1_unpaired.fastq.gz",
        r2_unpaired="results/trimmed/{sample}_R2_unpaired.fastq.gz"
    log:
        "results/logs/trimmomatic/{sample}.log"
    params:
        trimmer=config["trimmomatic"]["adapters"],
        extra=config["trimmomatic"]["extra_params"]
    threads: 4
    conda:
        "../envs/trimming.yaml"
    shell:
        """
        mkdir -p results/trimmed results/logs/trimmomatic
        trimmomatic PE -threads {threads} \
            {input.r1} {input.r2} \
            {output.r1} {output.r1_unpaired} \
            {output.r2} {output.r2_unpaired} \
            ILLUMINACLIP:{params.trimmer}:2:30:10 \
            {params.extra} 2> {log}
        """

rule cutadapt:
    """Alternative adapter trimming with cutadapt"""
    input:
        r1=lambda wildcards: f"data/{wildcards.sample}_R1{config['input']['fastq_suffix']}",
        r2=lambda wildcards: f"data/{wildcards.sample}_R2{config['input']['fastq_suffix']}"
    output:
        r1="results/cutadapt/{sample}_R1_trimmed.fastq.gz",
        r2="results/cutadapt/{sample}_R2_trimmed.fastq.gz"
    log:
        "results/logs/cutadapt/{sample}.log"
    params:
        adapters=config["cutadapt"]["adapters"],
        extra=config["cutadapt"]["extra_params"]
    threads: 4
    conda:
        "../envs/trimming.yaml"
    shell:
        "cutadapt -j {threads} {params.adapters} {params.extra} "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} > {log} 2>&1"