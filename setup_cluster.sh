# Example Snakemake profiles and configuration

# Profile for SLURM cluster execution
mkdir -p profiles/slurm

# Create cluster configuration
cat > profiles/slurm/config.yaml << 'EOF'
cluster: "sbatch --partition={resources.partition} --time={resources.time} --mem={resources.mem_mb} --cpus-per-task={threads} --job-name={rule} --output=logs/slurm/{rule}_%j.out --error=logs/slurm/{rule}_%j.err"
jobs: 100
use-conda: true
conda-frontend: mamba
latency-wait: 60
restart-times: 3
default-resources:
  - partition=general
  - time=4:00:00
  - mem_mb=8000
EOF

# Create resource specifications for rules
cat > profiles/slurm/resources.yaml << 'EOF'
bowtie2_build:
  mem_mb: 16000
  time: "2:00:00"

bowtie2_align:
  mem_mb: 12000
  time: "4:00:00"

macs2_callpeak:
  mem_mb: 8000
  time: "2:00:00"

bam_to_bigwig:
  mem_mb: 8000
  time: "1:00:00"
EOF