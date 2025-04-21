#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import pandas as pd


def parse_args():
    """
    Parse command-line arguments.

    --rIndex: Relative index for the BIM file.
    --phe: Trait identifier.
    """
    parser = argparse.ArgumentParser(description='Signal SNP annotation main script.')
    parser.add_argument('--rIndex', 
                        type=int, 
                        required=True,
                        help='Relative index in the BIM file to process.')
    
    parser.add_argument('--pheno', 
                        type=str, 
                        required=True,
                        help='Trait identifier (phenotype).')
    
    parser.add_argument(
        '--bfile',
        type=str,
        default="../example/small",
        help="Path to the plink genotype file."
    )
    
    parser.add_argument(
        '--HitFile',
        type=str,
        default="sig.index.summary.txt",
        help="Path storing the significant hits summary."
    )
    
    parser.add_argument(
        '--LDPath',
        type=str,
        default="fake_ld.bed",
        help="Path to LD information."
    )
    
    return parser.parse_args()


def main():
    # File paths for input data
    args = parse_args()
    rIndex = args.rIndex
    trait = args.pheno
    LDPath = args.LDPath
    plink_path = args.bfile
    HitFile = args.HitFile
    LD_filename = args.LDPath
    bim_filepath = f'{plink_path}.bim'

    # Load BIM file and assign column names
    bimfile = pd.read_csv(bim_filepath, sep='\t', header=None)
    bimfile.columns = ['chr', 'chrpos', 'MAF', 'pos', 'MAJ', 'MIN']

    # Initialize annotation arrays
    N = bimfile.shape[0]
    annot_local = np.zeros(N)
    annot_dist = np.zeros(N)
    annot_LD = np.zeros(N)
    data = []

    # Parse command-line arguments
    

    # Load BIM data again (alternative column names)
    df = pd.read_csv(bim_filepath, header=None, delimiter='\t')
    df.columns = ['chr', 'ID', 'unknown', 'loci', 'Ma', 'Mi']

    block = 100  # Block size (currently unused)

    # Load candidate index summary (kept for completeness)
    cand_df = pd.read_csv(HitFile,
                          header=None, delimiter=' ')

    # SNP information at the specified index
    inforow = bimfile.iloc[rIndex]
    CHR = inforow['chr']
    snpindex = inforow['pos']

    # Load LD region file for corresponding chromosome
    
    LD_pd = pd.read_csv(LD_filename, delim_whitespace=True)

    # Determine LD boundaries
    start_index = np.where(LD_pd.start <= snpindex)[0][-1]
    stop_index = np.where(LD_pd.stop >= snpindex)[0][0]
    start_pindex = LD_pd.start.iloc[start_index]
    stop_pindex = LD_pd.stop.iloc[stop_index]

    # Get locus and chromosome from df
    loci = df.iloc[rIndex, 3]
    Chr2 = df.iloc[rIndex, 0]

    # Construct sub-annotations
    subannot_local = annot_local.copy()
    subannot_local[df.chr == Chr2] = 1
    subannot_local[(df.chr == Chr2) & (df.loci >= start_pindex) & (df.loci <= stop_pindex)] = 0

    data.append([trait, rIndex, CHR, snpindex, start_pindex, stop_pindex])

    subannot_LD = annot_LD.copy()
    subannot_LD[(df.chr == Chr2) & (df.loci >= start_pindex) & (df.loci <= stop_pindex)] = 1

    subannot_dist = annot_dist.copy()
    subannot_dist[df.chr != Chr2] = 1

    subannot_rmLD = subannot_local + subannot_dist

    print(f'Index {rIndex}, local annotation ones: {np.sum(subannot_local)}; '
          f'distance annotation ones: {np.sum(subannot_dist)}')

    # Combine into final annotation matrix
    subannot = np.concatenate((
        subannot_LD.reshape(-1, 1),
        subannot_rmLD.reshape(-1, 1)
    ), axis=1)

    # Ensure output directory exists
    output_dir = 'annot'
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Save annotation and summary
    np.savetxt(f'{output_dir}/{trait}-{rIndex}.annot', subannot, fmt='%i')

    # LDsummary = pd.DataFrame(
    #     data=data,
    #     columns=['trait', 'Index', 'Chr', 'Pindex', 'Start_index', 'Stop_index']
    # )
    # LDsummary.to_csv(f'{output_dir}/LDsummary.csv', index=False)


if __name__ == "__main__":
    main()
