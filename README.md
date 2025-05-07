# FAME
This is the software for FAst Marginal Epistasis test (FAME) 

link: [https://www.biorxiv.org/content/10.1101/2023.09.10.557084v1.abstract]

## Prerequisites
The following packages are required on a Linux machine to compile and use the software package.
```
g++ (>=4.5)
cmake (>=2.8.12)
make (>=3.81)
```

## General description
Given a genotypes and phenotypes in PLINK format (bim, fam, bed, pheno), a target SNP (an index into the bim file), and an annotation file (with K annotations across the M SNPs), FAME jointly estimates the additive genetic variance across all the SNPs and the variance associated with the marginal epistatic effect (ME effect) of the target SNP. 

#### Note: 
Intuitively, the **annotation file** guides the software on how to partition the genotype matrix by SNPs. The column number of the annotation file determines the number of partitions. One basic rule of creating the annotation file: each row should contain one and only one "1"s.

To run the code, please delete the existing `build` folder (if any) and then use the following commands to compile the executable code.
```
mkdir build
cd build
cmake ..
make
```
The compilation time should be less than 1 minute. 


## Parameters
```
GENIE_GXG:
-g: path to genotype
-p: path to phenotype
-annot: annotation file where each SNP can be assigned to one of K annotations.
-gxgbin: the index of the annotation index to compute the ME effect. i.e. the set of SNPs that will be paired with the target SNP to define the marginal epistasis effect.  (These annotations are 0-based. So gxgbin of 0 refers to the SNPs with annotation of 1 in the first column in the annotation file. gxgbin of 1 refers to SNPs with annotation of 1 in the second column in the annotation file and so on.)
-snp: the index of the target SNP for which the marginal epistatic effect is to be estimated (This is a 1-based index that refers to the position of the SNP in the bim file)
-k: number of random vectors (we recommend 100 random vectors)
-jn: number of jackknife blocks
-o: output file
```

A full end-to-end pipeline to run FAME is provided in the `pipeline` directory. When running FAME, we recommend regressing out the effects of genetic variants that are in LD with the target SNP. This is handled by the pipeline. 
* To run the pipeline:
```
cd pipeline
bash run_fame_pipeline.sh
```

Two examples have been provided in the path `example`. 
* To test FAME with single bin:
```
cd example/single-bin
sh test.sh
```
* To test FAME with multiple bins:
```
cd example/multi-bin
sh test.sh
```

## Output:
FAME estimates K+2 variance components with K annotations corresponding to the annotation file. The first K components correspond to the additive component at each bin; component K+1 corresponds to the ME effect of the target SNP paired with the SNPs specified in gxgbin; and component K+2 corresponds to the i.i.d. noise component. 

Executing on the demo examples should take less than one minute to complete. 
