#!/bin/bash 
# set -x
## Input Arguments ###
genopath=${1:-../example/small}             # Path to PLINK genotype files (prefix)
phenopath=${2:-../example}                 # Path to input phenotype files
phenosnpfile=${3:-pheno-snp-pairs.txt}    # File containing pairs of phenotypes and SNP indices to be tested
code=${4:-../build/GENIE_GxG}              	# Executable binary for FAME
includecovar=false
covarflags=""

row_count=$(wc -l < ${phenosnpfile}) 		# Number of phenotype-index pairs

# Directory to store LD annotations
annotpath="annot"

# Directory to store residualized phenotypes
resphenopath="residualized_pheno"

# Directory to store FAME results
outpath="results"

# Number of random vectors for FAME
vec=100

gxgbin=1

if command -v python3 &>/dev/null; then
    PYTHON=python3
else
    PYTHON=python
fi

for phenosnppair in $(seq 1 $row_count); do
  
  echo "********** Running phenotype-SNP pair: $phenosnppair **********"
  
  var=$(awk "NR==${phenosnppair}" ${phenosnpfile})
  IFS=' ' read pheno snpid <<< $var # snpindex is the relative index in the plink file
  echo "SNP id: $snpid; Phenotype file: $pheno"

  # SNP index in bim file (0-index)
  snpindex=$(cat "${genopath}.bim" | awk -v a="$snpid" '$2==a { print NR - 1 }')
  if [ -z "$snpindex" ]; then
    echo "Cannot find $snpid in ${genopath}.bim"
    exit
  fi 
  echo "SNP index: $snpindex"

  # Original phenotype file
  phenofile="${phenopath}/${pheno}.pheno"

  if $includecovar ; then
    covarfile=`echo $phenofile|sed 's/.pheno$/.covar/'`
    if [ ! -f "$covarfile" ]; then
      echo "$covarfile does not exist"
      exit
    else
      echo "Using covariates in $covarfile"
	fi
	covarflags="-c $covarfile"
  fi

  # Path to annotation file 
  # containting LD block annotation relative to the target SNP
  # This is output by generate_ld_annotations.py
  annotfile="${annotpath}/${pheno}-${snpid}.annot"
	
  # Path to residualized phenotype
  # containing phenotype with genetic effects within the LD block residualized
  # This is output by linear_regression_annotation.py
  resphenofile="${resphenopath}/${pheno}-${snpid}.pheno"

  outfileprefix="${outpath}/${pheno}-${snpid}"
  outfile1="${outfileprefix}.res.out.txt"
  outfile2="${outfileprefix}.res.out.full.txt"

  if [ ! -d "${resphenopath}" ]; then
      echo "Creating directory: ${resphenopath}"
      mkdir -p "${resphenopath}"
  else
      echo "Directory ${resphenopath} exists"
  fi


  if [ ! -d "${annotpath}" ]; then
      echo "Creating directory: ${annotpath}"
      mkdir -p ${annotpath}
  else
      echo "Directory ${annotpath} exists"
  fi


  if [ ! -d "$outpath" ]; then
      echo "Creating directory: $outpath"
      mkdir -p $outpath
  else
      echo "Directory $outpath exists"
  fi

  if [ ! -f "$annotfile" ]; then 
    echo "$annotfile does not exist. Generating LD block annotations for target SNP ${snpindex} with respect to bim file ${genopath}.bim"
    $PYTHON "generate_ld_annotations.py" --rIndex ${snpindex} --output ${annotfile} 
  	echo "Generating residualized phenotype from ${phenofile}"
	$PYTHON "linear_regression_annotation.py" --rIndex ${snpindex} --bfile ${genopath} --annot ${annotfile} --pheno ${phenofile} --output ${resphenofile}
  	echo "Running FAME"
	$code  -g $genopath -p ${resphenofile} ${covarflags} -gxgbin ${gxgbin} -snp $((snpindex+1)) -k $vec  -jn 100  -o $outfile1 -annot $annotfile > $outfile2
  elif [ ! -f "$resphenofile" ]; then
    echo "$annotfile exists."
  	echo "Generating residualized phenotype from ${phenofile}"
	$PYTHON "linear_regression_annotation.py" --rIndex ${snpindex} --bfile ${genopath} --annot ${annotfile} --pheno ${phenofile} --output ${resphenofile}
  	echo "Running FAME"
	$code  -g $genopath -p ${resphenofile} ${covarflags} -gxgbin ${gxgbin} -snp $((snpindex+1)) -k $vec  -jn 100 -o $outfile1 -annot $annotfile > $outfile2
  elif [ ! -f "$outfile1" ]; then
  	echo "Running FAME"
	$code  -g $genopath -p ${resphenofile} ${covarflags} -gxgbin ${gxgbin} -snp $((snpindex+1)) -k $vec  -jn 100 -o $outfile1 -annot $annotfile > $outfile2
  else
  echo "Output file $outfile1 exists" 
  fi

done
