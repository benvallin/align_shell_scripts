#!/bin/bash

### Define help message
HELP='
NAME
  batch_pe_cutadapt.sh - run cutadapt iteratively on a set of PE samples

SYNOPSIS
  batch_pe_cutadapt.sh -s -o [-n -a -A -j -O -m -q -M -t]
  batch_pe_cutadapt.sh -h

DESCRIPTION
  Run cutadapt on a set of PE samples.
  -s  File listing the samples to operate on.
      Each line should contain the path to a sample-specific directory containing fastq files.
      Path can be absolute or relative to the current directory.
  -o  Path to output directory.
      Trimmed reads are written to sample-specific subdirectories within the provided directory.
      Path can be absolute or relative to the current directory.
  -n  find regexp matching the name of fastq files to trim. Default is \"*fastq.gz\".
      cutadapt expects paired-end files to be supplied in a pairwise fashion (e.g. file1_1.fq file1_2.fq file2_1.fq file2_2.fq).
      So, make sure the provided regex allows \"$(find path_to_sample_specific_directory -name regexp | sort)\" to return files in correct order.
  -a  Adapter sequence to be trimmed off read 1 (passed as -a to cutadapt).
  -A  Adapter sequence to be trimmed off read 2 (passed as -A to cutadapt).
  -j  The number of cores to be used for trimming (passed as -j to cutadapt). Default is 8.
  -O  Require minimum length overlap between read and adapter for an adapter to be found. (passed as -O to cutadapt). Default is 3.
  -m  Discard reads that became shorter than input value because of either quality or adapter trimming (passed as -m to cutadapt). Default is 0.
  -q  Quality cutoff score for trimming low-quality bases from 3 prime ends of each read before adapter removal. (passed as -q to cutadapt).
  -M  Discard reads with more Ns than provided number. Interpreted as a fraction of the read length if between 0 and 1. (passed as --max-n to cutadapt).
  -t  Trim Ns on ends of reads (passed as --trim-n to cutadapt).
  -h  Print this help.\n
'

### Parse option arguments
while getopts ":hs:o:n:a:A:j:O:m:q:M:t" opt
do
  case $opt in
    h) printf "$HELP"
    exit 0
    ;;
    s) SAMPLES="$OPTARG"
    ;;
    o) OUTDIR="${OPTARG%/}"
    ;;
    n) NAMES="$OPTARG"
    ;;
    a) ADAPTR1="$OPTARG"
    ;;
    A) ADAPTR2="$OPTARG"
    ;;
    j) THREADS="$OPTARG"
    ;;
    O) MINLENGTH="$OPTARG"
    ;;
    m) LENGTH="$OPTARG"
    ;;
    q) QUALITY="$OPTARG"
    ;;
    M) MAXN="$OPTARG"
    ;;
    t) TRIMN="--trim-n"
    ;;
    \?) printf "Invalid option -$OPTARG\n"
    exit 1
    ;;
    :) printf "Option -$OPTARG requires a valid argument\n" 
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) printf "Option -$opt requires a valid argument\n"
    exit 1
    ;;
  esac
done

# => Set names regexp to default value if not provided
if ! [[ -v NAMES ]]
then 
  NAMES="*fastq.gz"
fi

# => Set min length to default value if not provided
if ! [[ -v MINLENGTH ]]
then 
  MINLENGTH=3
fi

# => Set length to default value if not provided
if ! [[ -v LENGTH ]]
then 
  LENGTH=0
fi

# => Set n threads to default value if not provided
if ! [[ -v THREADS ]]
then 
  THREADS=8
fi

### Construct cutadapt's arguments
CUTADAPT_ARGS=()

# => Check that user provided adapters sequences and add them to cutadapt's arguments (-a and -A)
if [[ -v ADAPTR1 ]] && [[ -v ADAPTR2 ]]
then
  CUTADAPT_ARGS+=(-a $ADAPTR1 -A $ADAPTR2)
elif ! [[ -v ADAPTR1 ]] && ! [[ -v ADAPTR2 ]]
then
  printf "batch_pe_cutadapt needs user-provided adapters sequences. Please provide both -a and -A.\n"
  exit 1
elif ([[ -v ADAPTR1 ]] && ! [[ -v ADAPTR2 ]]) || (! [[ -v ADAPTR1 ]] && [[ -v ADAPTR2 ]])
then
  printf "batch_pe_cutadapt can only process PE reads. Please provide both -a and -A.\n"
  exit 1
fi

# => Add systematic arguments to cutadapt's arguments
CUTADAPT_ARGS+=(-O $MINLENGTH -m $LENGTH -j $THREADS)

# => Add -q VALUE, --max-n VALUE and --trim-n to cutadapt's arguments if provided 
if [[ -v QUALITY ]]
then
  CUTADAPT_ARGS+=(-q $QUALITY)
fi

if [[ -v MAXN ]]
then
  CUTADAPT_ARGS+=(--max-n $MAXN)
fi

if [[ -v TRIMN ]]
then
  CUTADAPT_ARGS+=($TRIMN)
fi

### Record samples details
N_SAMPLES=$(cat "$SAMPLES" | wc -l)
CURRENT_SAMPLE=0

printf "\n\n--->$N_SAMPLES SAMPLES TO PROCESS<---"

### Run cutadapt
for i in $(cat "$SAMPLES")
do 
  CURRENT_SAMPLE=$((CURRENT_SAMPLE+1))
  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES: PROCESSING SAMPLE $(basename $i)<---\n\n" 

  # => Create output directory if it does not exist already
  if ! [[ -d "$OUTDIR/$(basename $i)" ]]
  then 
    mkdir -p "$OUTDIR/$(basename $i)"
  fi

  cutadapt \
  "${CUTADAPT_ARGS[@]}" \
  -o "$OUTDIR"/$(basename "$i")/$(basename "$i")_R1.fastq.gz \
  -p "$OUTDIR"/$(basename "$i")/$(basename "$i")_R2.fastq.gz \
  $(find "$i" -name "$NAMES" | sort) 
    
  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES COMPLETE<---\n"
done

printf "\n\n--->ALL JOBS COMPLETE<---\n\n"


