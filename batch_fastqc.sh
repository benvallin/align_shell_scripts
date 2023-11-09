#!/bin/bash

### Define help message
HELP='
NAME
  batch_fastqc.sh - run fastqc on a set of samples

SYNOPSIS
  batch_fastqc.sh -s -o [-n -t]
  batch_fastqc.sh -h

DESCRIPTION
  Run fastqc on a set of samples.
  -s  File listing the samples to operate on.
      Each line should contain the path to a sample-specific directory containing fastq files.
      Path can be absolute or relative to the current directory.
  -o  Path to output directory (passed as -o to fastqc).
      Path can be absolute or relative to the current directory.
  -n  find regexp matching the name of fastq files to parse. Default is \"*fastq.gz\".
  -t  The number of threads to use concurrently (passed as -t to fastqc). Default is 8.
  -h  Print this help.\n
'

### Parse option arguments
while getopts ":hs:o:n:t:" opt
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
    t) THREADS="$OPTARG"
    ;;
    \?) printf "Invalid option -$OPTARG\n"
    exit 1
    ;;
    :) printf "Option -$OPTARG requires a valid argument\n" 
    exit 1
    ;;
  esac

  case $OPTARG in
    -*) printf "Option -$opt needs a valid argument\n"
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

### Record samples and fastq files details
N_SAMPLES=$(cat "$SAMPLES" | wc -l)

FQFILES=$(for i in $(cat "$SAMPLES"); do find "$i" -name "$NAMES"; done)

N_FQFILES=$(echo "$FQFILES" | wc -l)

printf "\n\n--->$N_FQFILES FASTQ FILES TO PROCESS OUT OF $N_SAMPLES SAMPLES<---\n\n"

### Create output directory if it does not exist already
if ! [[ -d "$OUTDIR" ]]
then 
  mkdir -p "$OUTDIR"
fi

### Run fastqc
fastqc -o "$OUTDIR" -t "$THREADS" $(echo "$FQFILES")

printf "\n\n--->ALL JOBS COMPLETE<---\n\n"

