"""
Peak Calling Rules
"""

def get_control_bam(wildcards):
    """Get control BAM file for a sample if it exists"""
    import pandas as pd
    samples_df = pd.read_csv(config["samples"], sep="\t")
    
    # Find the row for this sample
    sample_row = samples_df[samples_df['sample_id'] == wildcards.sample]
    if len(sample_row) == 0:
        return []
    
    control_id = sample_row['control'].iloc[0]
    
    # Check if control is not NaN/null and exists
    if pd.isna(control_id) or control_id == 'NA':
        return []
    else:
        return f"results/filtered/{control_id}.filtered.bam"

rule macs2_callpeak:
    """Call peaks using MACS2"""
    input:
        treatment="results/filtered/{sample}.filtered.bam",
        control=get_control_bam
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
        """
        if [ -n "{input.control}" ]; then
            macs2 callpeak -t {input.treatment} -c {input.control} -f BAMPE -g {params.genome_size} -n {params.name} --outdir {params.outdir} --broad {params.extra} > {log} 2>&1
        else
            macs2 callpeak -t {input.treatment} -f BAMPE -g {params.genome_size} -n {params.name} --outdir {params.outdir} --broad {params.extra} > {log} 2>&1
        fi
        """

def get_control_bedgraph(wildcards):
    """Get control bedgraph file for a sample if it exists"""
    import pandas as pd
    samples_df = pd.read_csv(config["samples"], sep="\t")
    
    # Find the row for this sample
    sample_row = samples_df[samples_df['sample_id'] == wildcards.sample]
    if len(sample_row) == 0:
        return []
    
    control_id = sample_row['control'].iloc[0]
    
    # Check if control is not NaN/null and exists
    if pd.isna(control_id) or control_id == 'NA':
        return []
    else:
        return f"results/bedgraph/{control_id}.bedgraph"

rule seacr_callpeak:
    """Call peaks using SEACR"""
    input:
        bedgraph="results/bedgraph/{sample}.bedgraph",
        control=get_control_bedgraph
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