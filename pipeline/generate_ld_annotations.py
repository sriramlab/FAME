#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import pandas as pd


def parse_args():
    """
    Parse command-line arguments.

    --rIndex: Relative index of target SNP in the BIM file.
    --output: Output file to save LD block annotations
    """
    parser = argparse.ArgumentParser(description='Target SNP annotation script.')
    parser.add_argument('--rIndex', 
                        type=int, 
                        required=True,
                        help='Relative index of the target SNP in the BIM file.')

    parser.add_argument(
        '--bfile',
        type=str,
        default="../example/small",
        help="Path to the plink genotype files (only need the bim file)"
    )
    
    parser.add_argument(
        '--LDPath',
        type=str,
        default="example_ld.bed",
        help="Path to LD information."
    )

    parser.add_argument(
        '--output',
        type=str,
        default="out.annot",
        help="Path to the output annotation file."
    )

    parser.add_argument(
        '--debug',
        action='store_true',
        help="Print debugging output"
    )


    return parser.parse_args()


def main():
    # File paths for input data
    args = parse_args()
    rIndex = args.rIndex
    LDPath = args.LDPath
    plink_path = args.bfile
    LD_filename = args.LDPath
    outfilename = args.output
    debug = args.debug
    bim_filepath = f'{plink_path}.bim'

    # Load BIM file and assign column names
    bimfile = pd.read_csv(bim_filepath, sep='\t', header=None)
    bimfile.columns = ['chr', 'chrpos', 'MAF', 'pos', 'MAJ', 'MIN']

    # Initialize annotation arrays
    N = bimfile.shape[0]
    annot_local = np.zeros(N)
    annot_dist = np.zeros(N)
    annot_LD = np.zeros(N)

    # Parse command-line arguments

    # SNP information at the specified index
    inforow = bimfile.iloc[rIndex]
    CHR = inforow['chr']
    snpindex = inforow['pos']
    print(f'Target SNP index {rIndex} corresponds to {CHR}:{snpindex} in {bim_filepath}')

    # Load LD region file for corresponding chromosome
    
    LD_pd = pd.read_csv(LD_filename, delim_whitespace=True) ## chr, block_start, block_end
    LD_pd['chr'] = LD_pd['chr'].astype(str)
    
    chr_mask = LD_pd['chr'] == str(CHR)
    assert snpindex <= LD_pd.loc[chr_mask,'stop'].max(), \
        f"Error: SNP index {snpindex} exceeds max stop index {LD_pd.stop.max()} in LD file."
    assert snpindex >= LD_pd.loc[chr_mask,'start'].min(), \
        f"Error: SNP index {snpindex} is less than min start index {LD_pd.start.min()} in LD file."
    
    # Determine LD boundaries
    start_index = np.where((LD_pd.chr==str(CHR))&(LD_pd.start <= snpindex))[0][-1]
    stop_index = np.where((LD_pd.chr==str(CHR))&(LD_pd.stop >= snpindex))[0][0]
    start_pindex = LD_pd.start.iloc[start_index]
    stop_pindex = LD_pd.stop.iloc[stop_index]
    
    assert start_pindex <= snpindex <= stop_pindex, \
        f"Error: SNP index {snpindex} not in range [{start_pindex}, {stop_pindex}]"
    

    # Get locus and chromosome from df
    loci = bimfile.iloc[rIndex, 3]
    Chr2 = bimfile.iloc[rIndex, 0]

    # Construct sub-annotations
    subannot_local = annot_local.copy()
    subannot_local[bimfile.chr == Chr2] = 1
    subannot_local[(bimfile.chr == Chr2) & (bimfile.pos >= start_pindex) & (bimfile.pos <= stop_pindex)] = 0

    subannot_LD = annot_LD.copy()
    subannot_LD[(bimfile.chr == Chr2) & (bimfile.pos >= start_pindex) & (bimfile.pos <= stop_pindex)] = 1

    subannot_dist = annot_dist.copy()
    subannot_dist[bimfile.chr != Chr2] = 1

    subannot_rmLD = subannot_local + subannot_dist

    if debug:
        print(f'Index {rIndex}, SNPs in LD block : {np.sum(subannot_LD)}; '
          f'SNPs outside LD block: {np.sum(subannot_rmLD)}')

    # Combine into final annotation matrix
    subannot = np.concatenate((
        subannot_LD.reshape(-1, 1),
        subannot_rmLD.reshape(-1, 1)
    ), axis=1)

    # The path for outfilename must exist. 
    # Currently path is created in run_fame_pipeline.sh
    np.savetxt(f'{outfilename}', subannot, fmt='%i')


if __name__ == "__main__":
    main()
