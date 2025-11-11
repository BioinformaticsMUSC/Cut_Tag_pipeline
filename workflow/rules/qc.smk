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