"""
Quality Control Rules
"""

rule fastqc:
    """Run FastQC on raw reads"""
    input:
        lambda wildcards: f"data/{wildcards.sample}_R{wildcards.read}{config['input']['fastq_suffix']}"
    output:
        html="results/fastqc/{sample}_R{read}_fastqc.html",
        zip="results/fastqc/{sample}_R{read}_fastqc.zip"
    params:
        outdir="results/fastqc"
    threads: 2
    conda:
        "../envs/qc.yaml"
    shell:
        """
        # Run FastQC
        fastqc -t {threads} -o {params.outdir} {input}
        
        # Rename outputs to match expected pattern (remove _001 if present)
        input_base=$(basename {input} {config[input][fastq_suffix]})
        if [[ "{config[input][fastq_suffix]}" == "_001.fastq.gz" ]]; then
            # Rename from sample_R1_001_fastqc.* to sample_R1_fastqc.*
            mv {params.outdir}/${{input_base}}_001_fastqc.html {output.html} 2>/dev/null || true
            mv {params.outdir}/${{input_base}}_001_fastqc.zip {output.zip} 2>/dev/null || true
        fi
        """

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

def get_multiqc_inputs():
    """Get MultiQC inputs based on whether trimming is skipped"""
    inputs = []
    # Always include raw FastQC
    inputs.extend(expand("results/fastqc/{sample}_R{read}_fastqc.zip", sample=SAMPLES, read=[1,2]))
    # Add trimmed FastQC only if trimming is not skipped
    if not config.get("processing", {}).get("skip_trimming", False):
        inputs.extend(expand("results/fastqc_trimmed/{sample}_R{read}_trimmed_fastqc.zip", sample=SAMPLES, read=[1,2]))
    # Always include BAM files
    inputs.extend(expand("results/aligned/{sample}_sorted.bam", sample=SAMPLES))
    return inputs

def get_multiqc_dirs():
    """Get MultiQC directories based on whether trimming is skipped"""
    dirs = ["results/fastqc", "results/aligned"]
    if not config.get("processing", {}).get("skip_trimming", False):
        dirs.append("results/fastqc_trimmed")
    return " ".join(dirs)

rule multiqc:
    """Generate MultiQC report"""
    input:
        get_multiqc_inputs()
    output:
        "results/qc/multiqc_report.html"
    params:
        outdir="results/qc",
        title="CUT&Tag Analysis Report",
        multiqc_dirs=get_multiqc_dirs()
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc -o {params.outdir} -n multiqc_report.html -i '{params.title}' {params.multiqc_dirs}"