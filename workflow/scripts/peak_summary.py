#!/usr/bin/env python3
"""
Simple peak summary without complex Bioconductor dependencies
"""

import pandas as pd
import os
from pathlib import Path

def summarize_peaks(peak_files, output_file):
    """
    Create a simple peak summary report
    """
    
    summary_data = []
    
    for peak_file in peak_files:
        sample_name = Path(peak_file).stem.replace('_peaks', '')
        
        # Count peaks
        with open(peak_file, 'r') as f:
            peak_count = sum(1 for line in f if not line.startswith('#') and line.strip())
        
        # Get basic stats if possible
        peaks_df = pd.read_csv(peak_file, sep='\t', header=None, comment='#')
        
        if len(peaks_df) > 0 and len(peaks_df.columns) >= 3:
            # Basic peak width statistics
            peak_widths = peaks_df.iloc[:, 2] - peaks_df.iloc[:, 1]
            
            summary_data.append({
                'sample': sample_name,
                'total_peaks': peak_count,
                'mean_width': peak_widths.mean() if len(peak_widths) > 0 else 0,
                'median_width': peak_widths.median() if len(peak_widths) > 0 else 0,
                'min_width': peak_widths.min() if len(peak_widths) > 0 else 0,
                'max_width': peak_widths.max() if len(peak_widths) > 0 else 0
            })
        else:
            summary_data.append({
                'sample': sample_name,
                'total_peaks': peak_count,
                'mean_width': 0,
                'median_width': 0,
                'min_width': 0,
                'max_width': 0
            })
    
    # Create summary dataframe
    summary_df = pd.DataFrame(summary_data)
    
    # Save summary
    summary_df.to_csv(output_file, sep='\t', index=False)
    
    print(f"Peak summary written to: {output_file}")
    print("\nSummary Statistics:")
    print(summary_df.to_string(index=False))

if __name__ == "__main__":
    peak_files = snakemake.input
    output_file = snakemake.output[0]
    
    summarize_peaks(peak_files, output_file)