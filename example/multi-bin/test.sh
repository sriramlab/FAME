#!/bin/bash

target='../target.txt'
gen=../small
phen=../h2_0.5.pheno
annot=multi-annot.txt
../../build/GENIE_GxG -g $gen -p $phen -gxgbin 1 -snp -1 -tp $target  -k 10 -jn 10  -o small.multi.results -annot $annot 






