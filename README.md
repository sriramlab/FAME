# FAME -- external target
This is the software extension for FAst Marginal Epistasis test (FAME) 

link: [https://www.biorxiv.org/content/10.1101/2023.09.10.557084v1.abstract]

## Prerequisites
The following packages are required on a Linux machine to compile and use the software package.
```
g++
cmakes
make
```

## General description
Given a genotype matrix (NxM) encoded as plink format (bim, fam, bed), the target trait in the same format, and an annotation file (MxK) with K bins across the M features, FAME jointly estimates the marginal effect across all the SNPs and the ME effect of the target SNP at k-th bin. 

#### Note: 
Intuitively, an **annotation file** guides the software on how to partition the genotype matrix by features. The column number of the annotation file determines the number of partitions made to the feature space. One basic rule of creating the annotation file: each row should contain one and only one "1", and all the rest columns of the same row should be 0s.

To run the code, please delete the existing `build` folder (if any) and then use the following commands to compile the executable code.
```
mkdir build
cd build
cmake ..
make
```



## Parameters
```
GENIE_GXG:
-g: path to genotype
-p: path to phenotype
-gxgbin: the bin index to compute the ME effect (0-based)
-snp: the target SNP index (1-based). If the target feature is provided outside of the genotype file, then put -1 and specify the target feature.
-tp: the external target feature to understand the ME effect. It should be a single column file, with N rows corresponding to each individual
-k: number of random vectors (100 is suggested)
-jn: number of jackknife blocks (100 is suggested)
-o: output file
-annot: annotation file
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

Compared to the default FAME version, the external target feature, whose marginal epistasis effect is of interest, needs to be provided in the standardized form (see `example/target.txt`) 

## Output:
K+2 estimated variance components will be reported. The first K components correspond to the additive component at each bin; the K+1 th component corresponds to the ME effect of the target bin; and the K+2 th component corresponds to the i.i.d. noise component. 
