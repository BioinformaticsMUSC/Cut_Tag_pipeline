#!/usr/bin/env python3
"""
Calculate TSS enrichment scores for CUT&Tag data quality assessment
"""

import pysam
import pandas as pd
import numpy as np
from pathlib import Path

def calculate_tss_enrichment(bam_file, tss_bed, window=2000):
    """
    Calculate TSS enrichment score
    TSS enrichment = (max signal in ±50bp around TSS) / (mean signal in flanking regions)
    """
    
    # Read TSS coordinates
    tss_sites = []
    with open(tss_bed, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3:
                chrom, start, end = fields[0], int(fields[1]), int(fields[2])
                tss_sites.append((chrom, (start + end) // 2))
    
    # Calculate coverage around TSS sites
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        tss_scores = []
        
        for chrom, tss_pos in tss_sites[:1000]:  # Limit to first 1000 TSS sites for speed
            try:
                # Get coverage in window around TSS
                start_pos = max(0, tss_pos - window)
                end_pos = tss_pos + window
                
                # Count reads in different regions
                center_reads = bam.count(chrom, tss_pos - 50, tss_pos + 50)  # ±50bp around TSS
                flank1_reads = bam.count(chrom, start_pos, tss_pos - 500)     # Upstream flank
                flank2_reads = bam.count(chrom, tss_pos + 500, end_pos)      # Downstream flank
                
                # Calculate enrichment
                center_density = center_reads / 100  # 100bp window
                flank_density = (flank1_reads + flank2_reads) / (2 * (window - 500))
                
                if flank_density > 0:
                    enrichment = center_density / flank_density
                    tss_scores.append(enrichment)
                    
            except Exception as e:
                continue
    
    if tss_scores:
        return np.mean(tss_scores)
    else:
        return 0.0

def main():
    bam_file = snakemake.input.bam
    tss_bed = snakemake.input.tss
    output_file = snakemake.output[0]
    
    sample_name = Path(bam_file).stem.replace('.filtered', '')
    
    enrichment_score = calculate_tss_enrichment(bam_file, tss_bed)
    
    # Write results
    with open(output_file, 'w') as f:
        f.write(f"sample\ttss_enrichment\n")
        f.write(f"{sample_name}\t{enrichment_score:.3f}\n")
    
    print(f"TSS enrichment score for {sample_name}: {enrichment_score:.3f}")

if __name__ == "__main__":
    main()