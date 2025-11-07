#!/usr/bin/env python3
"""
Generate consensus peaks across multiple samples
"""

import pandas as pd
from pathlib import Path
import subprocess
import tempfile

def merge_peaks(peak_files, output_file, min_overlap=0.5):
    """
    Merge peak files and create consensus peaks
    """
    
    if not peak_files:
        print("No peak files provided")
        return
    
    # Create temporary merged file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as temp_merged:
        temp_merged_path = temp_merged.name
    
    try:
        # Concatenate all peak files
        with open(temp_merged_path, 'w') as outf:
            for peak_file in peak_files:
                sample_name = Path(peak_file).stem.replace('_peaks', '')
                with open(peak_file, 'r') as inf:
                    for line in inf:
                        if not line.startswith('#') and line.strip():
                            fields = line.strip().split('\t')
                            if len(fields) >= 3:
                                # Add sample name as 4th column
                                outf.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\t{sample_name}\n")
        
        # Sort the merged file
        sorted_file = temp_merged_path + ".sorted"
        subprocess.run([
            "sort", "-k1,1", "-k2,2n", temp_merged_path
        ], stdout=open(sorted_file, 'w'), check=True)
        
        # Merge overlapping peaks
        merged_file = temp_merged_path + ".merged"
        subprocess.run([
            "bedtools", "merge", "-i", sorted_file, "-c", "4", "-o", "collapse"
        ], stdout=open(merged_file, 'w'), check=True)
        
        # Filter peaks present in at least 2 samples (or adjust threshold)
        min_samples = max(1, len(peak_files) // 2)  # At least half the samples
        
        with open(merged_file, 'r') as inf, open(output_file, 'w') as outf:
            outf.write("# Consensus peaks from CUT&Tag analysis\n")
            outf.write("# chrom\tstart\tend\tnum_samples\tsamples\n")
            
            for line in inf:
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    chrom, start, end, samples_str = fields[0], fields[1], fields[2], fields[3]
                    samples = samples_str.split(',')
                    unique_samples = list(set(samples))
                    
                    if len(unique_samples) >= min_samples:
                        outf.write(f"{chrom}\t{start}\t{end}\t{len(unique_samples)}\t{','.join(unique_samples)}\n")
        
        print(f"Generated consensus peaks in {output_file}")
        
        # Clean up temporary files
        Path(temp_merged_path).unlink(missing_ok=True)
        Path(sorted_file).unlink(missing_ok=True)
        Path(merged_file).unlink(missing_ok=True)
        
    except Exception as e:
        print(f"Error generating consensus peaks: {e}")
        # Clean up on error
        Path(temp_merged_path).unlink(missing_ok=True)

def main():
    peak_files = snakemake.input
    output_file = snakemake.output[0]
    
    merge_peaks(peak_files, output_file)

if __name__ == "__main__":
    main()