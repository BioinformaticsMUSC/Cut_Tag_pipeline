"""
Visualization Rules
"""

rule bam_to_bigwig:
    """Convert BAM files to BigWig for visualization"""
    input:
        bam="results/filtered/{sample}.filtered.bam",
        bai="results/filtered/{sample}.filtered.bam.bai"
    output:
        "results/bigwig/{sample}.bw"
    params:
        effective_genome_size=config["genome"]["effective_size"]
    threads: 4
    conda:
        "../envs/visualization.yaml"
    shell:
        "bamCoverage -b {input.bam} -o {output} "
        "--effectiveGenomeSize {params.effective_genome_size} "
        "--normalizeUsing RPKM "
        "--binSize 10 "
        "-p {threads}"

rule fragment_size_distribution:
    """Generate fragment size distribution plots"""
    input:
        expand("results/filtered/{sample}.filtered.bam", sample=SAMPLES)
    output:
        plot="results/qc/fragment_sizes.png",
        data="results/qc/fragment_sizes.txt"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/fragment_sizes.py"

rule tss_enrichment:
    """Calculate TSS enrichment scores"""
    input:
        bam="results/filtered/{sample}.filtered.bam",
        tss=config["genome"]["tss_bed"]
    output:
        "results/qc/{sample}_tss_enrichment.txt"
    conda:
        "../envs/visualization.yaml"
    script:
        "../scripts/tss_enrichment.py"

rule fingerprint_plot:
    """Generate BAM fingerprint plots for quality assessment"""
    input:
        expand("results/filtered/{sample}.filtered.bam", sample=SAMPLES)
    output:
        plot="results/qc/fingerprint_plot.png",
        metrics="results/qc/fingerprint_metrics.txt"
    params:
        labels=" ".join(SAMPLES)
    threads: 4
    conda:
        "../envs/visualization.yaml"
    shell:
        "plotFingerprint -b {input} "
        "--labels {params.labels} "
        "--minMappingQuality 30 "
        "--skipZeros "
        "--numberOfSamples 500000 "
        "-p {threads} "
        "--plotFile {output.plot} "
        "--outRawCounts {output.metrics}"

rule correlation_heatmap:
    """Generate sample correlation heatmap"""
    input:
        expand("results/bigwig/{sample}.bw", sample=SAMPLES)
    output:
        heatmap="results/qc/correlation_heatmap.png",
        matrix="results/qc/correlation_matrix.txt"
    params:
        labels=" ".join(SAMPLES)
    threads: 4
    conda:
        "../envs/visualization.yaml"
    shell:
        "multiBigwigSummary bins -b {input} "
        "--labels {params.labels} "
        "-o {output.matrix} "
        "-p {threads} && "
        "plotCorrelation -in {output.matrix} "
        "--corMethod pearson "
        "--skipZeros "
        "--plotTitle 'Pearson Correlation of Read Counts' "
        "--whatToPlot heatmap "
        "--colorMap RdYlBu "
        "--plotNumbers "
        "-o {output.heatmap}"