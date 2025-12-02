"""
Alignment Rules
"""

rule bowtie2_build:
    """Build Bowtie2 index for reference genome"""
    input:
        fasta=config["genome"]["fasta"]
    output:
        expand("resources/genome/genome_index.{ext}.bt2", ext=["1","2","3","4"]),
        "resources/genome/genome_index.rev.1.bt2",
        "resources/genome/genome_index.rev.2.bt2"
    params:
        prefix="resources/genome/genome_index"
    threads: 8
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2-build --threads {threads} {input.fasta} {params.prefix}"

def get_alignment_input_r1(wildcards):
    """Get R1 input file based on whether trimming is skipped"""
    if config.get("processing", {}).get("skip_trimming", False):
        return f"data/{wildcards.sample}_R1{config['input']['fastq_suffix']}"
    else:
        return f"results/trimmed/{wildcards.sample}_R1_trimmed.fastq.gz"

def get_alignment_input_r2(wildcards):
    """Get R2 input file based on whether trimming is skipped"""
    if config.get("processing", {}).get("skip_trimming", False):
        return f"data/{wildcards.sample}_R2{config['input']['fastq_suffix']}"
    else:
        return f"results/trimmed/{wildcards.sample}_R2_trimmed.fastq.gz"

rule bowtie2_align:
    """Align reads to reference genome (automatically chooses trimmed or raw reads)"""
    input:
        r1=get_alignment_input_r1,
        r2=get_alignment_input_r2,
        index=expand("resources/genome/genome_index.{ext}.bt2", ext=["1","2","3","4"])
    output:
        bam=temp("results/aligned/{sample}.bam"),
        log="results/logs/bowtie2/{sample}.log"
    params:
        index="resources/genome/genome_index",
        extra=config["bowtie2"]["extra_params"]
    threads: 8
    conda:
        "../envs/alignment.yaml"
    shell:
        "bowtie2 -x {params.index} -1 {input.r1} -2 {input.r2} "
        "-p {threads} {params.extra} 2> {output.log} | "
        "samtools view -Sb - > {output.bam}"

rule sort_bam:
    """Sort BAM files"""
    input:
        "results/aligned/{sample}.bam"
    output:
        "results/aligned/{sample}_sorted.bam"
    threads: 4
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools sort -@ {threads} -o {output} {input}"

rule index_bam:
    """Index sorted BAM files"""
    input:
        "results/aligned/{sample}_sorted.bam"
    output:
        "results/aligned/{sample}_sorted.bam.bai"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools index {input}"

rule alignment_stats:
    """Generate alignment statistics"""
    input:
        expand("results/aligned/{sample}_sorted.bam", sample=SAMPLES)
    output:
        "results/qc/alignment_stats.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        "echo 'Sample\tTotal_Reads\tMapped_Reads\tMapping_Rate' > {output}; "
        "for bam in {input}; do "
        "sample=$(basename $bam .sorted.bam); "
        "total=$(samtools view -c $bam); "
        "mapped=$(samtools view -c -F 4 $bam); "
        "rate=$(echo \"scale=2; $mapped/$total*100\" | bc); "
        "echo \"$sample\t$total\t$mapped\t$rate%\" >> {output}; "
        "done"