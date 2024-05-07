#!/bin/bash

target='../target.txt'
gen=../small
phen=../h2_0.5.pheno
annot=single.annot.small.txt
../../build/GENIE_GxG -g $gen -p $phen -gxgbin 0 -snp -1 -tp $target  -k 10 -jn 10  -o small.single.results -annot $annot 






