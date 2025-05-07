#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import pandas as pd
from bed_reader import open_bed
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm

def parse_args():
    # Parse command-line arguments.
    parser = argparse.ArgumentParser(
        description="Process trait configuration, trait ID, and annotation file."
    )

    parser.add_argument(
        '--annot',
        type=str,
        required=True,
        help="Path to the annotation file containing SNP annotation data."
    )
    
    
    parser.add_argument(
        '--bfile',
        type=str,
        default="../example/small",
        help="Path to the plink genotype file."
    )
    
    parser.add_argument(
        '--pheno',
        type=str,
        default="../example/h2_0.5.pheno",
        help="Path to the ivrted phenotype file."
    )

    parser.add_argument(
        '--output',
        type=str,
        default="out.pheno",
        help="Path to the residualized output phenotype file."
    )

    parser.add_argument(
        '--debug',
        action='store_true',
        help="Print debugging output"
    )

    return parser.parse_args()
    

def main():
    
    args = parse_args()
    # Set variables based on command-line arguments.
    annot_path = args.annot
    bfile = args.bfile
    base_phe_path = args.pheno
    debug = args.debug
    outfilename = args.output
	
    # Load annotation data.
    annot = np.loadtxt(annot_path)
    
    trait_config = base_phe_path.split('/')[-1]

    # Define file paths and load genotype data.
    bimfile = pd.read_csv(f'{bfile}.bim', delim_whitespace=True, header=None)
    bimfile.columns = ['CHR', 'BP', 'x', 'SNP', 'y', 'z']
    
    bed_file = f'{bfile}.bed'
    G = open_bed(bed_file)

    # Load phenotype data.
    base_phe_df = pd.read_csv(base_phe_path, delimiter=' ')

    if debug:
        print(base_phe_df.head())

    base_phe = base_phe_df.pheno.values

    # Logging message.
    if debug:
       print(f'######### In trait {trait_config} #########')

    # Identify indices where the annotation is equal to 1.
    ld_indices = np.where(annot[:, 0] == 1)[0]
    if debug:
       print(f'LD indices: {ld_indices}')
    istart = ld_indices[0]
    iend = ld_indices[-1]

    # Read genotype data for the LD region.
    XLD = G.read(np.s_[:, istart:(iend + 1)])
    imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
    XLD = imputer.fit_transform(XLD)
    
    # Remove genetic component from phenotype.
    newPhe_res = base_phe - LinearRegression().fit(XLD, base_phe).predict(XLD)

    # Determine target SNP index from the annotation file name.
    annot_filename = os.path.basename(annot_path)
    target_index = int(annot_filename.split('-')[-1].split('.')[0])
    if debug:
       print(f'target index: {target_index}')

    # Read and process the target SNP's genotype data.
    X_t = G.read(np.s_[:, target_index:(target_index + 1)])
    imputer = SimpleImputer(missing_values=np.nan, strategy='mean')
    X_t = imputer.fit_transform(X_t.reshape(-1, 1)).flatten()
    X_t = sm.add_constant(X_t)

    # Fit the linear regression model using statsmodels.
    model = sm.OLS(newPhe_res, X_t)
    result = model.fit()
    coefficients = result.params
    standard_errors = result.bse
    print(f'p-value is {result.pvalues}')

    # Save the updated phenotype results.
    newPhe_res_df = base_phe_df.copy()
    newPhe_res_df['pheno'] = newPhe_res
    newPhe_res_df.iloc[:, 0] = newPhe_res_df.iloc[:, 0].astype(int) ## FID
    newPhe_res_df.iloc[:, 1] = newPhe_res_df.iloc[:, 1].astype(int) ## IID


    # The path for outfilename must exist. 
    # Currently path is created in run_fame_pipeline.sh
    newPhe_res_df.to_csv(f'{outfilename}', sep=' ', index=None)
    print (f'Writing residualized phenotypes to {outfilename}')

    if debug:
       print(f'annot : {annot_filename}')
       print(f'outfilename : {outfilename}')

if __name__ == "__main__":
    main()

