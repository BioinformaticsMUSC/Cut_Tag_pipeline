"""
Peak Calling Rules
"""

rule macs2_callpeak:
    """Call peaks using MACS2"""
    input:
        treatment="results/filtered/{sample}.filtered.bam",
        control=lambda wildcards: f"results/filtered/{samples_df[samples_df['sample_id'] == wildcards.sample]['control'].iloc[0]}.filtered.bam" if not samples_df[samples_df['sample_id'] == wildcards.sample]['control'].isna().iloc[0] else []
    output:
        narrowpeak="results/peaks/macs2/{sample}_peaks.narrowPeak",
        summits="results/peaks/macs2/{sample}_summits.bed",
        broadpeak="results/peaks/macs2/{sample}_peaks.broadPeak",
        xls="results/peaks/macs2/{sample}_peaks.xls"
    params:
        name="{sample}",
        outdir="results/peaks/macs2",
        genome_size=config["genome"]["effective_size"],
        extra=config["macs2"]["extra_params"]
    log:
        "results/logs/macs2/{sample}.log"
    conda:
        "../envs/peaks.yaml"
    shell:
        "macs2 callpeak "
        "-t {input.treatment} "
        + ("-c {input.control} " if len(input.control) > 0 else "") +
        "-f BAMPE "
        "-g {params.genome_size} "
        "-n {params.name} "
        "--outdir {params.outdir} "
        "--broad "
        "{params.extra} > {log} 2>&1"

rule seacr_callpeak:
    """Call peaks using SEACR"""
    input:
        bedgraph="results/bedgraph/{sample}.bedgraph",
        control=lambda wildcards: f"results/bedgraph/{samples_df[samples_df['sample_id'] == wildcards.sample]['control'].iloc[0]}.bedgraph" if not samples_df[samples_df['sample_id'] == wildcards.sample]['control'].isna().iloc[0] else []
    output:
        stringent="results/peaks/seacr/{sample}.stringent.bed",
        relaxed="results/peaks/seacr/{sample}.relaxed.bed"
    params:
        prefix="results/peaks/seacr/{sample}",
        threshold=config["seacr"]["threshold"]
    conda:
        "../envs/peaks.yaml"
    script:
        "../scripts/run_seacr.py"

rule bam_to_bedgraph:
    """Convert BAM to bedgraph for SEACR"""
    input:
        bam="results/filtered/{sample}.filtered.bam",
        bai="results/filtered/{sample}.filtered.bam.bai"
    output:
        "results/bedgraph/{sample}.bedgraph"
    conda:
        "../envs/peaks.yaml"
    shell:
        "bedtools genomecov -ibam {input.bam} -bg > {output}"

rule peak_annotation:
    """Annotate peaks with genomic features"""
    input:
        peaks="results/peaks/macs2/{sample}_peaks.narrowPeak",
        gtf=config["genome"]["gtf"]
    output:
        "results/peaks/annotated/{sample}_annotated_peaks.txt"
    conda:
        "../envs/peaks.yaml"
    script:
        "../scripts/annotate_peaks.R"

rule peak_consensus:
    """Generate consensus peaks across samples"""
    input:
        expand("results/peaks/macs2/{sample}_peaks.narrowPeak", sample=SAMPLES)
    output:
        "results/peaks/consensus_peaks.bed"
    conda:
        "../envs/peaks.yaml"
    script:
        "../scripts/consensus_peaks.py"