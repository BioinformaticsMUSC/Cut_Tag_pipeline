#!/usr/bin/env python3
"""
Run SEACR peak calling
"""

import subprocess
import sys
from pathlib import Path

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
    
    # Run SEACR in stringent mode
    if control_bedgraph:
        cmd_stringent = [
            "bash", str(seacr_script),
            treatment_bedgraph, control_bedgraph,
            "non", "stringent", f"{output_prefix}.stringent"
        ]
        cmd_relaxed = [
            "bash", str(seacr_script), 
            treatment_bedgraph, control_bedgraph,
            "non", "relaxed", f"{output_prefix}.relaxed"
        ]
    else:
        cmd_stringent = [
            "bash", str(seacr_script),
            treatment_bedgraph, str(threshold),
            "non", "stringent", f"{output_prefix}.stringent"
        ]
        cmd_relaxed = [
            "bash", str(seacr_script),
            treatment_bedgraph, str(threshold), 
            "non", "relaxed", f"{output_prefix}.relaxed"
        ]
    
    print(f"Running SEACR stringent: {' '.join(cmd_stringent)}")
    subprocess.run(cmd_stringent, check=True)
    
    print(f"Running SEACR relaxed: {' '.join(cmd_relaxed)}")
    subprocess.run(cmd_relaxed, check=True)

if __name__ == "__main__":
    treatment = snakemake.input.bedgraph
    output_prefix = snakemake.params.prefix
    control = snakemake.input.control[0] if len(snakemake.input.control) > 0 else None
    threshold = snakemake.params.threshold
    
    run_seacr(treatment, output_prefix, control, threshold)