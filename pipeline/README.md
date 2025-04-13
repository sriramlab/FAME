# FAME Main Pipeline

This bash script is a demo **end-to-end pipeline** for:

1. **Generating SNP-level LD annotations**
2. **Residualizing phenotypes using local LD SNPs**
3. **Running FAME GxG interaction analysis** using a compiled binary (`GENIE_GxG`)

The pipeline processes a list of significant SNPs defined in a summary file and performs **per-SNP annotation, regression, and model evaluation**.

---

## Pipeline Overview

For each SNP (`Index`) and associated phenotype (`pheno`), `estimation_res.sh` did the following:

1. **Generate annotation file** from LD structure:
    - via `5.split_LD.py`
2. **Residualize phenotype** by removing genetic variance from the LD block:
    - via `5.5.reg_Lin.py`
3. **Run GxG model** to evaluate interactions using the residual phenotype:
    - via `GENIE_GxG` binary

To run the pipeline, simply executing:
```
bash estimation_res.sh
```

Results will be generated at
```
results/xxx.pheno.res.out.txt ## core results
results/xxx.pheno.res.out.full.txt ## Full results
```
The example output is
```
sigma^2_0: -0.0281539 se: 0.0036717 ## var_g of LD block -- can be safely discarded
sigma^2_1: 0.454532 se: 0.133152 ## var_g of non-LD block
sigma^2_2: -0.0066988 se: 0.0526437 ## var_gxg of non-LD block
sigma^2_3: 0.526866 se: 0.0945873 ## var_e

```

# Generate annotation file
`5.split.LD.py`
This script generates SNP-level annotations based on local linkage disequilibrium (LD) structure.


Given:
- A PLINK `.bim` genotype file,
- A significant hit index (i.e., a SNP relative index in the bim file),
- A phenotype label,
- An LD region definition file (fake LD info file provided, the original LD block info can be downloaded [here](https://bitbucket.org/nygcresearch/ldetect-data/src/master/)),
- An optional list of significant hits,

The script produces:
- A two-column `.annot` file containing (0-based):
  - Column 1: LD region indicator (1 if SNP is in the same LD block as the hit)
  - Column 2: Non-LD region where LD will be computed upon (1 if SNP is on a different chromosome or outside the local LD block)


# Residualize phenotype
`5.5.reg_Lin.py`
This script processes genotype, phenotype, and annotation data to **residualize phenotype values** by removing the genetic component explained by SNPs in an LD block. It then fits a univariate linear regression at a **target SNP** to estimate its association with the residualized phenotype.

Given:
- A **PLINK genotype** file (`.bed`, `.bim`),
- An **annotation file** marking LD regions and a target SNP,
- A **phenotype file** (PLINK-style format),

The script:
1. Reads SNPs in an LD region defined by the annotation file.
2. Removes the genetic signal of those SNPs from the phenotype using linear regression.
4. Outputs the updated residualized phenotype in PLINK format.
