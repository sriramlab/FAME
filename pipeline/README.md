# FAME Main Pipeline

This bash script is a demo **end-to-end pipeline** to run FAME. This pipeline performs the following steps :

1. **Generates SNP-level annotations to indicate the LD block around the target SNP**
2. **Residualizes the phenotype w.r.t additive genetic effects of SNPs in the LD block**
3. **Runs FAME GxG interaction analysis** using a compiled binary (`GENIE_GxG`)

The pipeline processes a list of target SNP-trait pairs and performs **annotation of SNPs that lie in the LD block of the target SNP, regression of additive genetic effects within this block, and marginal espistasis analysis of each target SNP-trait pair**.

---

## Pipeline Overview

To test for ME at a target SNP (indicated by its index in a PLINK bim file) on a phenotype, run `run_fame_pipeline.sh` to perform each of these steps:

1. **Generate SNP-level annotation file** from LD blocks:
    - via `generate_ld_annotations.py`
2. **Residualize phenotype** by regressing the additive effects of genetic variants that lie within the LD block:
    - via `linear_regression_annotation.py`
3. **Run FAME** to test for marginal espistasis using the residual phenotype:
    - via `FAME/build/GENIE_GxG` binary

Inputs
```
genopath=../example/small    ##Path to PLINK genotype files (prefix)
phenopath=../example/ ## Path to input phenotype files 
phenosnpfile=pheno-snp-pairs.txt ## File containing pairs of phenotypes and SNP indices to be tested
code=../build/GENIE_GxG  ## Executable binary for FAME
# Note:
# Phenotype values are expected to be properly normalized (e.g., via inverse-rank normal transformation, IVRT) before running the pipeline.
# If the executable ../build/GENIE_GxG is run separately, the -snp argument expects a one-based index.
```

To run the pipeline, run:
```
bash run_fame_pipeline.sh
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
`generate_ld_annotations.py`
This script generates SNP-level annotations based on linkage disequilibrium (LD) blocks.

Given:
- A target SNP (specified by an index in the bim file), 
- A path to a PLINK `.bim` genotype file (example PLINK files provided in ../example/small.{bed,bim,fam}), 
- A bed file defining the LD blocks (example LD info file provided in example_ld.bed, the original LD block info can be downloaded [here](https://bitbucket.org/nygcresearch/ldetect-data/src/master/)),

The script produces:
- A two-column `.annot` file where:
  - Column 1 is 1 if SNP is in the same LD block as the target SNP
  - Column 2 is 1 if SNP is on a different chromosome or outside the LD block containing the target SNP

IMPORTANT: The index of the target SNP is 0-based.

# Residualize phenotype
`linear_regression_annotation.py`
This script processes genotype, phenotype, and annotation data to **residualize phenotype values** by removing the genetic component explained by SNPs in the LD block corresponding to a target SNP. As a check, it then fits a univariate linear regression at the **target SNP** to estimate its association with the residualized phenotype (this step should produce a non-significant association).

Given:
- A **PLINK genotype** file (`.bed`, `.bim`),
- An **annotation file** marking the LD block containing the target SNP,
- A **phenotype file** (PLINK-style format),

The script:
1. Reads SNPs in an LD region defined by the annotation file.
2. Regresses out those SNPs from the phenotype using linear regression.
3. Outputs the residualized phenotype in PLINK format (.pheno).
