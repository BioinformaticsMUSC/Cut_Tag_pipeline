"""
Peak Calling Rules
"""

def get_control_bam(wildcards):
    """Get control BAM file for a sample if it exists"""
    import csv
    
    with open(config["samples"], 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['sample_id'] == wildcards.sample:
                control_id = row['control']
                # Check if control is not empty/NA and exists
                if control_id and control_id != 'NA':
                    return f"results/filtered/{control_id}.filtered.bam"
                else:
                    return []
    return []

def is_control_sample(sample_id):
    """Check if a sample is a control (has NA in control column)"""
    import csv
    
    with open(config["samples"], 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['sample_id'] == sample_id:
                control_id = row['control']
                return not control_id or control_id == 'NA'
    return False

def get_treatment_samples():
    """Get list of samples that are not controls"""
    import csv
    treatment_samples = []
    
    with open(config["samples"], 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            control_id = row['control']
            # Include sample if it has a control (not NA or empty)
            if control_id and control_id != 'NA':
                treatment_samples.append(row['sample_id'])
    
    return treatment_samples

rule macs2_callpeak:
    """Call peaks using MACS2 (skip control samples)"""
    input:
        treatment="results/filtered/{sample}.filtered.bam",
        control=get_control_bam
    output:
        narrowpeak="results/peaks/macs2/{sample}_peaks.narrowPeak",
        summits="results/peaks/macs2/{sample}_summits.bed",
        #broadpeak="results/peaks/macs2/{sample}_peaks.broadPeak",
        xls="results/peaks/macs2/{sample}_peaks.xls"
    params:
        name="{sample}",
        outdir="results/peaks/macs2",
        genome_size=config["genome"]["effective_size"],
        extra=config["macs2"]["extra_params"],
        samples_file=config["samples"]
    log:
        "results/logs/macs2/{sample}.log"
    conda:
        "../envs/peaks.yaml"
    shell:
        """
        # Check if this is a control sample using awk/grep
        # Look for the sample in the TSV file and check if control column is NA
        CONTROL_VALUE=$(awk -F'\t' -v sample="{wildcards.sample}" '$1 == sample {{print $4}}' {params.samples_file} | head -1)
        
        if [ "$CONTROL_VALUE" = "NA" ] || [ -z "$CONTROL_VALUE" ]; then
            # Create empty output files for control samples
            mkdir -p {params.outdir}
            touch {output.narrowpeak} {output.summits} {output.xls}
            echo "Skipped peak calling for control sample {wildcards.sample}" > {log}
            echo "Skipped peak calling for control sample {wildcards.sample}"
        else
            # Run peak calling for treatment samples
            if [ -n "{input.control}" ]; then
                macs2 callpeak -t {input.treatment} -c {input.control} -f BAMPE -g {params.genome_size} -n {params.name} --outdir {params.outdir} {params.extra} > {log} 2>&1
            else
                macs2 callpeak -t {input.treatment} -f BAMPE -g {params.genome_size} -n {params.name} --outdir {params.outdir} {params.extra} > {log} 2>&1
            fi
        fi
        """

def get_control_bedgraph(wildcards):
    """Get control bedgraph file for a sample if it exists"""
    import csv
    
    with open(config["samples"], 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['sample_id'] == wildcards.sample:
                control_id = row['control']
                # Check if control is not empty/NA and exists
                if control_id and control_id != 'NA':
                    return f"results/bedgraph/{control_id}.bedgraph"
                else:
                    return []
    return []

rule seacr_callpeak:
    """Call peaks using SEACR (skip control samples)"""
    input:
        bedgraph="results/bedgraph/{sample}.bedgraph",
        control=get_control_bedgraph
    output:
        stringent="results/peaks/seacr/{sample}.stringent.bed",
        relaxed="results/peaks/seacr/{sample}.relaxed.bed"
    params:
        prefix="results/peaks/seacr/{sample}",
        threshold=config["seacr"]["threshold"],
        samples_file=config["samples"]
    conda:
        "../envs/peaks.yaml"
    script:
        "../scripts/run_seacr_conditional.py"

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

# Commented out due to complex Bioconductor dependencies
# Uncomment and install annotation packages separately if needed
# rule peak_annotation:
#     """Annotate peaks with genomic features"""
#     input:
#         peaks="results/peaks/macs2/{sample}_peaks.narrowPeak",
#         gtf=config["genome"]["gtf"]
#     output:
#         "results/peaks/annotated/{sample}_annotated_peaks.txt"
#     conda:
#         "../envs/peaks_annotation.yaml"  # Would need separate environment
#     script:
#         "../scripts/annotate_peaks.R"

rule peak_summary:
    """Generate simple peak summary statistics (treatment samples only)"""
    input:
        lambda wildcards: expand("results/peaks/macs2/{sample}_peaks.narrowPeak", sample=get_treatment_samples())
    output:
        "results/peaks/peak_summary.txt"
    conda:
        "../envs/peaks.yaml"
    script:
        "../scripts/peak_summary.py"

rule peak_consensus:
    """Generate consensus peaks across treatment samples"""
    input:
        lambda wildcards: expand("results/peaks/macs2/{sample}_peaks.narrowPeak", sample=get_treatment_samples())
    output:
        "results/peaks/consensus_peaks.bed"
    conda:
        "../envs/peaks.yaml"
    script:
        "../scripts/consensus_peaks.py"