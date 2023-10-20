# FAME
This is the software for FAst Marginal Epistasis test (FAME)

Given a genotype matrix (NxM) encoded as plink format (bim, fam, bed), the target trait in the same format, and an annotation file (MxK) with K bins across the M features, FAME jointly estimates the marginal effect across all the SNPs and the ME effect of the target SNP at k-th bin. 

To run the code, please delete the existing build (if any) and then use the following commands to compile the executable code.
```
mkdir build
cd build
cmake ..
make
```

Commands
```
GENIE_GXG:
-g: path to genotype
-p: path to phenotype
-gxgbin: the bin index to compute the ME effect (0-based)
-snp: the target SNP index (1-based)
-k: number of random vector
-jn: number of jackknife blocks
-o: output file
--annot: annotation file
```

