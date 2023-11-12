#!/bin/bash

### Define help message
HELP='
NAME
  batch_pe_cogent_demux.sh - run cogent demux iteratively on a set of PE samples 

SYNOPSIS
  batch_pe_cogent_demux.sh -s -b -o [-t -1 -2 -n -u -S -g]
  batch_pe_cogent_demux.sh -h

DESCRIPTION
  Run cogent demux on a set of PE samples.
  -s  File listing the samples to operate on.
      Each line should contain the path to a sample-specific directory containing pair-end reads fastq files.
      Path can be absolute or relative to the current directory.
  -b  Path to well list file (passed as -b to cogent demux). 
      Path can be absolute or relative to the current directory.
  -o  Path to output directory.
      Results are written to barcode-specific subdirectories within the provided directory.
      Path can be absolute or relative to the current directory.
  -t  Format string describing the experimental protocol used (passed as -t to cogent demux). Default is \"ICELL8_FLA\".
  -1  find regexp matching the name of fastq files containing the #1 mates. Default is \"*_R1_*\".
  -2  find regexp matching the name of fastq files containing the #1 mates. Default is \"*_R2_*\".
  -n  The number of demultiplexing processes to use concurrently (passed as -n to cogent demux). Default is 8.
  -u  Save undetermined/unselected/short reads to undetermined FASTQ files (passed as --undetermined_fq to cogent demux).
  -S  Output merged FASTQ files (passed as --no_split_fastqs to cogent demux). 
      Barcodes are written into read names and merged into large FASTQ file. 
  -g  Do not compress (gzip) output FASTQ files (passed as --no_gz to cogent demux).
  -h  Print this help.\n
'

### Parse option arguments
while getopts ":hs:b:o:t:1:2:n:uSg" opt
do
  case $opt in
    h) printf "$HELP"
    exit 0
    ;;
    s) SAMPLES="$OPTARG"
    ;;
    b) WELLS="$OPTARG"
    ;;
    o) OUTDIR="${OPTARG%/}"
    ;;
    t) EXPTYPE="$OPTARG"
    ;;
    1) MATES1="$OPTARG"
    ;;
    2) MATES2="$OPTARG"
    ;;
    n) THREADS="$OPTARG"
    ;;
    u) UNDETERMINED="--undetermined_fq"
    ;;
    S) NOSPLIT="--no_split_fastqs"
    ;;
    g) NOGZ="--no_gz"
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

# => Set experimental protocol type to default value if not provided
if ! [[ -v EXPTYPE ]]
then 
  EXPTYPE="ICELL8_FLA"
fi

# => Set MATES1 regexp to default value if not provided
if ! [[ -v MATES1 ]]
then 
  MATES1="*_R1_*"
fi

# => Set MATES2 regexp to default value if not provided
if ! [[ -v MATES2 ]]
then 
  MATES2="*_R2_*"
fi

# => Set n threads to default value if not provided
if ! [[ -v THREADS ]]
then 
  THREADS=8
fi

### Construct Cogent's arguments
COGENT_ARGS=()

# => Add systematic arguments to Cogent's arguments
COGENT_ARGS+=(-b $WELLS -t $EXPTYPE -n $THREADS)

# => Add --undetermined_fq, --no_split_fastqs and --no_gz to Cogent's arguments if provided 
OPTIONAL_ARGS="$UNDETERMINED $NOSPLIT $NOGZ"

for i in "$OPTIONAL_ARGS"
do COGENT_ARGS+=($i)
done

### Record samples details
N_SAMPLES=$(cat "$SAMPLES" | wc -l)
CURRENT_SAMPLE=0

printf "\n\n--->$N_SAMPLES SAMPLES TO PROCESS<---"

### Run cogent demux
for i in $(cat "$SAMPLES")
do 
  CURRENT_SAMPLE=$((CURRENT_SAMPLE+1))
  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES: PROCESSING SAMPLE $(basename $i)<---\n\n"
  
  cogent demux "${COGENT_ARGS[@]}" \
  -i $(find "$i" -name "$MATES1") \
  -p $(find "$i" -name "$MATES2") \
  -o "$OUTDIR"/$(basename "$i")/

  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES COMPLETE<---\n"
done

printf "\n\n--->ALL JOBS COMPLETE<---\n\n"
