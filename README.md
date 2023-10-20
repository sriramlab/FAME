# FAME
This is the software for FAst Marginal Epistasis test (FAME) 

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

To run the code, please delete the existing build (if any) and then use the following commands to compile the executable code.
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
-snp: the target SNP index (1-based)
-k: number of random vector
-jn: number of jackknife blocks
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

## Output:
K+2 estimated variance components will be reported. The first K components correspond to the additive component at each bin; the K+1 th component corresponds to the ME effect of the target bin; and the K+2 th component corresponds to the i.i.d. noise component.
