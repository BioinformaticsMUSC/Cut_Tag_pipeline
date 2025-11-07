#!/usr/bin/env python3
"""
Generate fragment size distribution plots for CUT&Tag data
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import pandas as pd
import numpy as np
from pathlib import Path

def get_fragment_sizes(bam_file, max_reads=100000):
    """Extract fragment sizes from paired-end BAM file"""
    fragments = []
    
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        count = 0
        for read in bam.fetch():
            if count >= max_reads:
                break
            
            # Only count properly paired reads with positive insert size
            if (read.is_proper_pair and 
                read.template_length > 0 and 
                read.template_length < 1000):  # Filter out very large fragments
                fragments.append(read.template_length)
                count += 1
    
    return fragments

def plot_fragment_distribution(bam_files, output_plot, output_data):
    """Plot fragment size distributions for all samples"""
    
    all_fragments = {}
    
    # Collect fragment sizes from all BAM files
    for bam_file in bam_files:
        sample_name = Path(bam_file).stem.replace('.filtered', '')
        print(f"Processing {sample_name}...")
        
        fragments = get_fragment_sizes(bam_file)
        all_fragments[sample_name] = fragments
    
    # Create combined dataframe
    df_list = []
    for sample, fragments in all_fragments.items():
        df_list.append(pd.DataFrame({
            'fragment_size': fragments,
            'sample': sample
        }))
    
    df = pd.concat(df_list, ignore_index=True)
    
    # Save raw data
    df.to_csv(output_data, sep='\t', index=False)
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Plot 1: Histogram overlay
    for sample in all_fragments.keys():
        sample_data = df[df['sample'] == sample]['fragment_size']
        ax1.hist(sample_data, bins=50, alpha=0.7, label=sample, density=True)
    
    ax1.set_xlabel('Fragment Size (bp)')
    ax1.set_ylabel('Density')
    ax1.set_title('Fragment Size Distribution')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Box plot
    sns.boxplot(data=df, x='sample', y='fragment_size', ax=ax2)
    ax2.set_xlabel('Sample')
    ax2.set_ylabel('Fragment Size (bp)')
    ax2.set_title('Fragment Size Distribution by Sample')
    ax2.tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Print summary statistics
    print("\nFragment Size Summary:")
    print(df.groupby('sample')['fragment_size'].describe())

if __name__ == "__main__":
    bam_files = snakemake.input
    output_plot = snakemake.output.plot
    output_data = snakemake.output.data
    
    plot_fragment_distribution(bam_files, output_plot, output_data)