# CUT&Tag Analysis Pipeline

Snakemake pipeline for CUT&Tag (Cleavage Under Targets and Tagmentation) data analysis.

## Overview

This pipeline processes CUT&Tag sequencing data through the following steps:

1. **Quality Control**: FastQC analysis of raw and trimmed reads
2. **Read Trimming**: Adapter and quality trimming with Trimmomatic/Cutadapt
3. **Alignment**: Mapping to reference genome using Bowtie2
4. **Filtering**: Quality filtering and duplicate removal
5. **Peak Calling**: Peak identification using MACS2 and SEACR
6. **Annotation**: Peak annotation with genomic features
7. **Visualization**: BigWig generation and quality plots
8. **Quality Assessment**: Fragment size analysis, TSS enrichment, correlation analysis

## Directory Structure

```
cut_tag/
├── workflow/
│   ├── Snakefile              # Main workflow
│   ├── rules/                 # Individual rule modules
│   │   ├── qc.smk
│   │   ├── trimming.smk
│   │   ├── alignment.smk
│   │   ├── filtering.smk
│   │   ├── peaks.smk
│   │   └── visualization.smk
│   ├── scripts/               # Custom scripts
│   └── envs/                  # Conda environments
├── config/
│   ├── config.yaml           # Pipeline configuration
│   └── samples.tsv           # Sample information
├── data/                     # Input FASTQ files
├── resources/                # Reference genome and annotations
├── results/                  # Output files
└── README.md                 # This file
```

## Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/) >= 7.0
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)
- Sufficient disk space for reference genomes and results

## Installation

1. Clone this repository:
```bash
git clone <repository-url>
cd cut_tag
```

2. Install Snakemake (if not already installed):
```bash
conda install -c conda-forge -c bioconda snakemake
```

## Quick Start

### 1. Prepare Input Data

Place your paired-end FASTQ files in the `data/` directory. The pipeline supports flexible naming conventions:

**Default pattern (with `_001` suffix):**
```
{sample_id}_R1_001.fastq.gz
{sample_id}_R2_001.fastq.gz
```

**Alternative patterns:** Modify `fastq_suffix` in `config/config.yaml`:
- For files without `_001`: change to `.fastq.gz`
- For other patterns: adjust accordingly

### 2. Configure Samples

Edit `config/samples.tsv` with your sample information:
```tsv
sample_id	condition	replicate	control	antibody
H3K4me3_rep1	H3K4me3	1	IgG_rep1	H3K4me3
H3K27ac_rep1	H3K27ac	1	IgG_rep1	H3K27ac
IgG_rep1	IgG	1	NA	IgG
```

### 3. Configure Pipeline

Edit `config/config.yaml` to specify:
- Reference genome files
- Analysis parameters
- Quality thresholds

### 4. Prepare Reference Genome

Download and prepare your reference genome files:
```bash
# Example for mouse mm10
mkdir -p resources/genome
cd resources/genome

# Download genome fasta
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
mv mm10.fa genome.fa

# Download GTF annotation
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz
gunzip gencode.vM25.annotation.gtf.gz
mv gencode.vM25.annotation.gtf genome.gtf
```

### 5. Run the Pipeline

Run a dry-run to check the workflow:
```bash
snakemake --dry-run --cores 1
```

Execute the pipeline:
```bash
# Local execution
snakemake --cores 8 --use-conda

# With resource limits
snakemake --cores 8 --use-conda --resources mem_mb=32000

# On a cluster (example for SLURM)
snakemake --profile slurm --cores 100
```

## Configuration Options

### Key Parameters in `config.yaml`:

**Reference Genome:**
- `genome.fasta`: Reference genome FASTA file
- `genome.gtf`: Gene annotation GTF file  
- `genome.effective_size`: Effective genome size for MACS2

**Quality Filtering:**
- `filtering.min_quality`: Minimum mapping quality (default: 30)
- `filtering.max_fragment_size`: Maximum fragment size for nucleosome-free reads (default: 120)

**Peak Calling:**
- `macs2.extra_params`: Additional MACS2 parameters
- `seacr.threshold`: FDR threshold for SEACR peak calling

## Output Files

### Key Results:

- `results/qc/multiqc_report.html`: Comprehensive QC report
- `results/peaks/macs2/{sample}_peaks.narrowPeak`: MACS2 peaks
- `results/peaks/seacr/{sample}.stringent.bed`: SEACR peaks
- `results/bigwig/{sample}.bw`: Signal tracks for visualization
- `results/peaks/consensus_peaks.bed`: Consensus peaks across samples

### Directory Structure:
```
results/
├── fastqc/                   # FastQC reports
├── trimmed/                  # Trimmed FASTQ files
├── aligned/                  # BAM files
├── filtered/                 # Quality-filtered BAMs
├── peaks/                    # Peak calling results
│   ├── macs2/
│   ├── seacr/
│   └── annotated/
├── bigwig/                   # BigWig signal tracks
└── qc/                       # Quality control reports
```

## Quality Control Metrics

The pipeline generates several QC metrics:

1. **Fragment Size Distribution**: Nucleosome-free vs nucleosome-bound fragments
2. **TSS Enrichment**: Signal enrichment around transcription start sites
3. **Library Complexity**: Assessment of PCR duplication rates
4. **Peak Statistics**: Number and distribution of identified peaks
5. **Sample Correlation**: Correlation between replicates

## Advanced Usage

### Custom Peak Calling

To use only MACS2 or SEACR, modify the `all` rule in `workflow/Snakefile`.

### Adding New Samples

1. Add FASTQ files to `data/`
2. Update `config/samples.tsv`
3. Re-run the pipeline (Snakemake will process only new samples)

### Parallel Execution

For large datasets, use cluster execution:
```bash
# Create cluster profile
mkdir -p ~/.config/snakemake/slurm
# Add cluster configuration files

# Submit to cluster
snakemake --profile slurm --jobs 100
```

## Troubleshooting

### Common Issues:

1. **Memory errors**: Increase `--resources mem_mb=` parameter
2. **Missing conda packages**: Update environment YAML files
3. **Reference genome errors**: Verify file paths in config.yaml
4. **Peak calling failures**: Check input BAM file quality and parameters

### Logs and Debugging:

- Check `results/logs/` for detailed error messages
- Use `snakemake --reason` to understand rule execution
- Run with `--verbose` for detailed output

## Citation

If you use this pipeline, please cite:

- **Snakemake**: Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2012.
- **CUT&Tag Protocol**: Kaya-Okur et al. "CUT&Tag for efficient epigenomic profiling of small samples and single cells". Nature Communications 2019.

## License

This pipeline is released under the MIT License.

## Support

For questions and issues:
1. Check the troubleshooting section above
2. Review Snakemake documentation
3. Open an issue in this repository