"""
BAM Filtering Rules
"""

rule filter_bam:
    """Filter BAM files for high-quality reads"""
    input:
        bam="results/aligned/{sample}.sorted.bam",
        bai="results/aligned/{sample}.sorted.bam.bai"
    output:
        "results/filtered/{sample}.filtered.bam"
    params:
        quality=config["filtering"]["min_quality"],
        flags=config["filtering"]["sam_flags"]
    threads: 4
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools view -@ {threads} -b -q {params.quality} {params.flags} "
        "{input.bam} > {output}"

rule remove_duplicates:
    """Remove PCR duplicates using Picard"""
    input:
        "results/filtered/{sample}.filtered.bam"
    output:
        bam="results/filtered/{sample}.dedup.bam",
        metrics="results/qc/picard/{sample}.dedup_metrics.txt"
    conda:
        "../envs/alignment.yaml"
    shell:
        "picard MarkDuplicates "
        "INPUT={input} "
        "OUTPUT={output.bam} "
        "METRICS_FILE={output.metrics} "
        "REMOVE_DUPLICATES=true "
        "VALIDATION_STRINGENCY=LENIENT"

rule filter_fragments:
    """Filter fragments by size (keep nucleosome-free reads < 120bp)"""
    input:
        "results/filtered/{sample}.filtered.bam"
    output:
        "results/filtered/{sample}.nuc_free.bam"
    params:
        max_size=config["filtering"]["max_fragment_size"]
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools view -h {input} | "
        "awk 'BEGIN{{OFS=\"\\t\"}} "
        "$1~/^@/ || ($9>0 && $9<{params.max_size}) || ($9<0 && $9>-{params.max_size})' | "
        "samtools view -b - > {output}"

rule index_filtered_bam:
    """Index filtered BAM files"""
    input:
        "results/filtered/{sample}.filtered.bam"
    output:
        "results/filtered/{sample}.filtered.bam.bai"
    conda:
        "../envs/alignment.yaml"
    shell:
        "samtools index {input}"