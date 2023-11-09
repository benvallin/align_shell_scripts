#!/bin/bash

### Define help message
HELP='
NAME
  dataset_to_sample_directory.sh - move sample-specific fastq files from input directory to sample-specific directories  

SYNOPSIS
  dataset_to_sample_directory.sh -i -o [-n -p]
  dataset_to_sample_directory.sh -h

DESCRIPTION
  Move sample-specific fastq files stored in input directory to sample-specific directories.
  -i  Path to input directory containing the fastq files to move. Can be absolute or relative to the current directory. 
      File names should contain the ID of the corresponding sample. 
  -o  Path to output directory. Can be absolute or relative to the current directory.
      Files are moved to sample-specific subdirectories within the provided directory.
  -n  find regexp matching the name of fastq files to move. Default is \"*fastq.gz\".
  -p  sed regexp matching a string to remove from the file names in order to generate sample-specific names. Default is \"_R[12].*$\".
  -h  Print this help.\n
'

### Parse option arguments
while getopts ":hi:o:n:p:" opt
do
  case $opt in
    h) printf "$HELP"
    exit 0
    ;;
    i) INDIR="${OPTARG%/}"
    ;;
    o) OUTDIR="${OPTARG%/}"
    ;;
    n) NAMES="$OPTARG"
    ;;
    p) PATTERN="$OPTARG"
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

# => Set pattern regexp to default value if not provided
if ! [[ -v PATTERN ]]
then 
  PATTERN="_R[12].*$"
fi

### Record samples details
FILES=$(for i in $(find "$INDIR" -name "$NAMES" | sort); do basename "$i"; done)

SAMPLES=$(sed "s/$PATTERN//g" <(echo "$FILES") | sort | uniq)

N_SAMPLES=$(echo "$SAMPLES" | wc -l)

printf "\n\n--->$N_SAMPLES SAMPLES TO PROCESS<---\n"

### Move fastq files to sample-specific directories
for i in $(echo "$SAMPLES")
do 
  FILESI=$(for j in $(echo "$FILES" | grep -e "$i"); do echo "$INDIR/$j"; done)
  printf "\n\nSample $i:\nMoving:\n$FILESI\nTo:\n$OUTDIR/$i/\n"
  if ! [[ -d "$OUTDIR"/"$i" ]]
  then 
    mkdir -p "$OUTDIR"/"$i"
  fi
  mv $(echo "$FILESI") "$OUTDIR"/"$i"
done

printf "\n--->JOB COMPLETE<---\n\n"
