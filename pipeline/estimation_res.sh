#!/bin/bash

for SGE_TASK_ID in $(seq 1 2); do
  
  echo "SGE_TASK_ID: $SGE_TASK_ID"

  gen=../example/small ## genotype plink file

  code=../build/GENIE_GxG ## FAME code


  dircp="traits/${phe}_res"
  out=results/
  vec=100

  sumfile="sig.index.summary.txt"
  var=$(awk "NR==${SGE_TASK_ID}" ${sumfile})
  IFS=' ' read pheno Index <<< $var # Index is the relative index in the plink file
  echo "index is $Index; pheno: $pheno"

  annotfile=annot/${pheno}-${Index}.annot
  gxgbin=1


  dircp="traits/${pheno}_res"
  echo "create directory $dircp"


  if [ ! -d "results" ]; then
      mkdir -p "results"
  fi
  echo "executing file $dircp/${pheno}-${Index}.annot"

  if [ ! -f "$annotfile" ]; then 
    echo "$annotfile not exists -- start generating"
    python "5.split_LD.py" --pheno ${pheno} --rIndex ${Index} --HitFile ${sumfile} &&
    python "5.5.reg_Lin.py" --phe_path ../example/${pheno} --bfile ${gen} --annot ${annotfile} &&
    $code  -g $gen -p $dircp/${pheno}-${Index}.annot  -gxgbin ${gxgbin} -snp $((Index+1)) -k $vec  -jn 100  -o $out/${Index}.${pheno}.res.out.txt  -annot $annotfile > $out/${Index}.${pheno}.res.out.full.txt
    sleep 2
    rm "$dircp/${pheno}-${Index}.annot"
  elif [ ! -f "$dircp/${pheno}-${Index}.annot" ]; then
      echo "generate res phe"
      python "5.5.reg_Lin.py" --phe_path ../example/${pheno} --bfile ${gen} --annot ${annotfile} &&
      $code  -g $gen -p $dircp/${pheno}-${Index}.annot  -gxgbin ${gxgbin} -snp $((Index+1)) -k $vec  -jn 100  -o $out/${Index}.${pheno}.res.out.txt  -annot $annotfile > $out/${Index}.${pheno}.res.out.full.txt
      sleep 2
      rm "$dircp/${pheno}-${Index}.annot"
  elif [ ! -f "$out/${Index}.${pheno}.res.out.txt" ]; then
    $code  -g $gen -p $dircp/${pheno}-${Index}.annot  -gxgbin ${gxgbin} -snp $((Index+1)) -k $vec  -jn 100  -o $out/${Index}.${pheno}.res.out.txt  -annot $annotfile > $out/${Index}.${pheno}.res.out.full.txt
    sleep 2
    rm "$dircp/${pheno}-${Index}.annot"
  else
    echo "$out/${Index}.${pheno}.res.out.txt exists" 
  fi

done