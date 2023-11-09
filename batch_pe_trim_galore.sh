#!/bin/bash

### Define help message
HELP='
NAME
  batch_pe_trim_galore.sh - run trim_galore iteratively on a set of PE samples

SYNOPSIS
  batch_pe_trim_galore.sh -s -o [-n -a -A -I -N -S -j -q -O -m -M -t]
  batch_pe_trim_galore.sh -h

DESCRIPTION
  Run trim_galore on a set of PE samples.
  -s  File listing the samples to operate on.
      Each line should contain the path to a sample-specific directory containing fastq files.
      Path can be absolute or relative to the current directory.
  -o  Path to output directory.
      Trimmed reads are written to sample-specific subdirectories within the provided directory.
      Path can be absolute or relative to the current directory.
  -n  find regexp matching the name of fastq files to trim. Default is \"*fastq.gz\".
      trim_galore expects paired-end files to be supplied in a pairwise fashion (e.g. file1_1.fq file1_2.fq file2_1.fq file2_2.fq).
      So, make sure the provided regex allows \"$(find path_to_sample_specific_directory -name regexp | sort)\" to return files in correct order.
  -a  Adapter sequence to be trimmed off read 1 (passed as -a to trim_galore).
      If not specified, Trim Galore will try to auto-detect the adapter type.
      See option \"-a\" in \"trim_galore --help\". Equivalent to option \"-a\" in \"cutadapt --help\". 
  -A  Adapter sequence to be trimmed off read 2 (passed as -a2 to trim_galore).
      If not specified, Trim Galore will try to auto-detect the adapter type.
      See option \"-a2\" in \"trim_galore --help\". Equivalent to option \"-A\" in \"cutadapt --help\". 
  -I  Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter 'AGATCGGAAGAGC' (passed as --illumina to trim_galore).
      See option \"--illumina\" in \"trim_galore --help\".
  -N  Adapter sequence to be trimmed is the first 12bp of the Nextera adapter 'CTGTCTCTTATA' (passed as --nextera to trim_galore).
      See option \"--nextera\" in \"trim_galore --help\".
  -S  Adapter sequence to be trimmed off read 1 is the first 12bp of the Illumina Small RNA 3p adapter 'TGGAATTCTCGG'.
      Adapter sequence to be trimmed off read 2 is the first 12bp of the Illumina Small RNA 5p adapter 'GATCGTCGGACT'.
      Selecting to trim Small RNA adapters will also lower the -m value to 18bp (passed as --small_rna to trim_galore). 
      See option \"--small_rna\" in \"trim_galore --help\".
  -j  The number of cores to be used for trimming (passed as -j to trim_galore). Equivalent to option \"-j\" in \"cutadapt --help\". Default is 8.
  -q  Trim low-quality ends from reads in addition to adapter removal (passed as -q to trim_galore).
      See option \"-q\" in \"trim_galore --help\". Relates to (but may differ from) option \"-q\" in \"cutadapt --help\"). Default is 20.
  -O  Overlap with adapter sequence required to trim a sequence (passed as --stringency to trim_galore).
      See option \"--stringency\" in \"trim_galore --help\". Equivalent to option \"-O\" in \"cutadapt --help\". Default is 3.
  -m  Discard reads that became shorter than input value because of either quality or adapter trimming (passed as --length to trim_galore).
      See option \"--length\" in \"trim_galore --help\". Equivalent to option \"-m\" in \"cutadapt --help\". Default is 20.
  -M  The total number of Ns a read may contain before it will be removed altogether (passed as --max_n to trim_galore).
      See option \"--max_n\" in \"trim_galore --help\". Equivalent to option \"--max-n\" in \"cutadapt --help\". 
  -t  Removes Ns from either side of the read (passed as --trim-n to trim_galore).
      See option \"--trim-n\" in \"trim_galore --help\". Equivalent to option \"--trim-n\" in \"cutadapt --help\".
  -h  Print this help.\n
'

### Parse option arguments
while getopts ":hs:o:n:a:A:INSj:q:O:m:M:t" opt
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
    I) ILLUMINA="--illumina"
    ;;
    N) NEXTERA="--nextera"
    ;;
    S) SMALL_RNA="--small_rna"
    ;;
    j) THREADS="$OPTARG"
    ;;
    q) QUALITY="$OPTARG"
    ;;
    O) STRINGENCY="$OPTARG"
    ;;
    m) LENGTH="$OPTARG"
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

# => Set n threads to default value if not provided
if ! [[ -v THREADS ]]
then 
  THREADS=8
fi

# => Set quality to default value if not provided
if ! [[ -v QUALITY ]]
then 
  QUALITY=20
fi

# => Set stringency to default value if not provided
if ! [[ -v STRINGENCY ]]
then 
  STRINGENCY=3
fi

# => Set length to default value if not provided
if ! [[ -v LENGTH ]]
then 
  LENGTH=20
fi

### Construct trim_galore's arguments
TRIM_GALORE_ARGS=()

# => Determine if:
# - user provided actual adapters sequences (-a and -A)
# - user provided an adapters type option (--illumina or --nextera or --small_rna)
# - adapters type should be auto-detected (no option)
# => Add user-provided options to trim_galore's arguments when necessary
if [[ -v ADAPTR1 ]] && [[ -v ADAPTR2 ]]
then
  ADAPT_MODE="SEQUENCES PROVIDED BY USER"
  TRIM_GALORE_ARGS+=(-a $ADAPTR1 -a2 $ADAPTR2)
elif ([[ -v ADAPTR1 ]] && ! [[ -v ADAPTR2 ]]) || (! [[ -v ADAPTR1 ]] && [[ -v ADAPTR2 ]])
then
  printf "batch_pe_trim_galore can only process PE reads. Please provide both -a and -A.\n"
  exit 1
fi

TOTAL_INS=0
for i in ILLUMINA NEXTERA SMALL_RNA 
do 
  if [[ -v "$i" ]]; then TOTAL_INS=$((TOTAL_INS+1)); fi
done

if [[ $ADAPT_MODE == "SEQUENCES PROVIDED BY USER" ]] && ! [[ $TOTAL_INS -eq 0 ]]
then
  printf "Please provide either adapters sequences or one of \"-I\", \"-N\" or \"-S\", not both.\n"
  exit 1
elif ! [[ $ADAPT_MODE == "SEQUENCES PROVIDED BY USER" ]] && [[ $TOTAL_INS -gt 1 ]]
then
  printf "Please provide only one of \"-I\", \"-N\" or \"-S\", not several.\n"
  exit 1
elif ! [[ $ADAPT_MODE == "SEQUENCES PROVIDED BY USER" ]] && [[ $TOTAL_INS -eq 1 ]]
then
  ADAPT_MODE="OPTION \"$ILLUMINA$NEXTERA$SMALL_RNA\" PROVIDED BY USER"
  TRIM_GALORE_ARGS+=($ILLUMINA$NEXTERA$SMALL_RNA)
elif ! [[ $ADAPT_MODE == "SEQUENCES PROVIDED BY USER" ]] && [[ $TOTAL_INS -eq 0 ]]
then
  ADAPT_MODE="AUTO-DETECTION"
fi

# => Add systematic arguments to trim_galore's arguments
TRIM_GALORE_ARGS+=(--paired --retain_unpaired --quality $QUALITY --stringency $STRINGENCY --length $LENGTH --cores $THREADS)

# => Add --max_n VALUE and --trim-n to trim_galore's arguments if provided 
if [[ -v MAXN ]]
then
  TRIM_GALORE_ARGS+=(--max_n $MAXN)
fi

if [[ -v TRIMN ]]
then
  TRIM_GALORE_ARGS+=($TRIMN)
fi

### Record samples details
N_SAMPLES=$(cat "$SAMPLES" | wc -l)
CURRENT_SAMPLE=0

printf "\n\n--->$N_SAMPLES SAMPLES TO PROCESS<---"

printf "\n\n--->ADAPTER MODE: $ADAPT_MODE<---"

### Run trim_galore
for i in $(cat "$SAMPLES")
do 
  CURRENT_SAMPLE=$((CURRENT_SAMPLE+1))
  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES: PROCESSING SAMPLE $(basename $i)<---\n\n" 

  trim_galore "${TRIM_GALORE_ARGS[@]}" -o "$OUTDIR"/$(basename "$i")/ $(find "$i" -name "$NAMES" | sort)
    
  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES COMPLETE<---\n"
done

printf "\n\n--->ALL JOBS COMPLETE<---\n\n"


