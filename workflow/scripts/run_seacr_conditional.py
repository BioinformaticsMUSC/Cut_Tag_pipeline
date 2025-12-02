#!/usr/bin/env python3
"""
Run SEACR peak calling with conditional logic for control samples
"""

import subprocess
import sys
import csv
from pathlib import Path

def is_control_sample(sample_id, samples_file):
    """Check if a sample is a control (has NA in control column)"""
    with open(samples_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['sample_id'] == sample_id:
                control_id = row['control']
                return not control_id or control_id == 'NA'
    return False

def run_seacr(treatment_bedgraph, output_prefix, control_bedgraph=None, threshold=0.01):
    """
    Run SEACR peak calling
    """
    
    # Download SEACR if not available
    seacr_script = Path("workflow/scripts/SEACR_1.3.sh")
    if not seacr_script.exists():
        print("Downloading SEACR...")
        subprocess.run([
            "wget", "-O", str(seacr_script),
            "https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.sh"
        ], check=True)
        subprocess.run(["chmod", "+x", str(seacr_script)], check=True)
    
    seacr_R_script = Path("workflow/scripts/SEACR_1.3.R")
    if not seacr_R_script.exists():
        print("Downloading SEACR R script...")
        subprocess.run([
            "wget", "-O", str(seacr_R_script),
            "https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.R"
        ], check=True)
        subprocess.run(["chmod", "+x", str(seacr_R_script)], check=True)
        
    # Run SEACR in stringent mode
    if control_bedgraph:
        cmd_stringent = [
            "bash", str(seacr_script),
            treatment_bedgraph, control_bedgraph,
            "non", "stringent", f"{output_prefix}"
        ]
        cmd_relaxed = [
            "bash", str(seacr_script), 
            treatment_bedgraph, control_bedgraph,
            "non", "relaxed", f"{output_prefix}"
        ]
    else:
        cmd_stringent = [
            "bash", str(seacr_script),
            treatment_bedgraph, str(threshold),
            "non", "stringent", f"{output_prefix}"
        ]
        cmd_relaxed = [
            "bash", str(seacr_script),
            treatment_bedgraph, str(threshold), 
            "non", "relaxed", f"{output_prefix}"
        ]
    
    print(f"Running SEACR stringent: {' '.join(cmd_stringent)}")
    subprocess.run(cmd_stringent, check=True)
    
    print(f"Running SEACR relaxed: {' '.join(cmd_relaxed)}")
    subprocess.run(cmd_relaxed, check=True)

if __name__ == "__main__":
    # Get sample ID from output prefix
    output_prefix = snakemake.params.prefix
    sample_id = Path(output_prefix).name
    
    # Check if this is a control sample
    samples_file = snakemake.params.samples_file
    
    if is_control_sample(sample_id, samples_file):
        # Create empty output files for control samples
        print(f"Skipping SEACR peak calling for control sample {sample_id}")
        Path(snakemake.output.stringent).touch()
        Path(snakemake.output.relaxed).touch()
    else:
        # Run SEACR for treatment samples
        treatment = snakemake.input.bedgraph
        output_prefix = snakemake.params.prefix
        
        # Handle control input - could be empty list or string
        control = None
        if hasattr(snakemake.input, 'control') and snakemake.input.control:
            if isinstance(snakemake.input.control, list):
                control = snakemake.input.control[0] if len(snakemake.input.control) > 0 else None
            else:
                control = snakemake.input.control
        
        threshold = snakemake.params.threshold
        
        run_seacr(treatment, output_prefix, control, threshold)