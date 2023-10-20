# FAME
This is the software for FAst Marginal Epistasis test (FAME)

Given a genotype matrix (NxM) encoded as plink format (bim, fam, bed), the target trait in the same format, and an annotation file (MxK) with K bins across the M features, FAME jointly estimate the marginal effect across all the SNPs and the ME effect of the target SNP at k-th bin. 

To run the code, please first use the following commands.
```
mkdir build
cd build
cmake ..
make
```
