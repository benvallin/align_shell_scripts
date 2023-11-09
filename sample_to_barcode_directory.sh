#!/bin/bash

### Define help message
HELP='
NAME
  sample_to_barcode_directory.sh - move barcode-specific fastq files from sample-specific to barcode-specific directories  

SYNOPSIS
  sample_to_barcode_directory.sh -s -o [-p -d]
  sample_to_barcode_directory.sh -h

DESCRIPTION
  Move barcode-specific fastq files stored in sample-specific directories to barcode-specific directories.
  -s  File listing the samples to operate on.
      Each line should contain the path to a sample-specific directory containing barcode-specific fastq files.
      Path can be absolute or relative to the current directory. File names should contain the ID of the corresponding barcode. 
  -o  Path to output directory. Can be absolute or relative to the current directory.
      Files are moved to barcode-specific subdirectories within the provided directory.
  -p  sed regexp matching a string to remove from the file names in order to generate barcode-specific names. Default is \"_R[12].*$\".
  -d  If provided, the input sample-specific directories are deleted once the job is complete.
  -h  Print this help.\n
'

### Parse option arguments
while getopts ":hs:o:p:d" opt
do
  case $opt in
    h) printf "$HELP"
    exit 0
    ;;
    s) SAMPLES="$OPTARG"
    ;;
    o) OUTDIR="${OPTARG%/}"
    ;;
    p) PATTERN="$OPTARG"
    ;;
    d) DELETE=
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

# => Set pattern regexp to default value if not provided
if ! [[ -v PATTERN ]]
then 
  PATTERN="_R[12].*$"
fi

### Record samples details
N_SAMPLES=$(cat "$SAMPLES" | wc -l)

printf "\n\n--->$N_SAMPLES SAMPLES TO PROCESS<---\n"

### Move fastq files to barcode-specific directories
for i in $(cat "$SAMPLES")
do 
  BARCODES=$(sed "s/$PATTERN//g" <(ls "$i"))
  N_BARCODES=$(echo "$BARCODES" | wc -l)
  printf "\n--->Sample $(basename $i): creating $N_BARCODES barcode-specific directories<---" 
  for j in $(echo "$BARCODES")
  do
    if ! [[ -d "$OUTDIR"/"$j" ]]
    then 
      mkdir -p "$OUTDIR"/"$j"
    fi
    for k in $(find "$i" -name "$j*")
    do 
      mv "$k" "$OUTDIR"/"$j"
    done
  done
done

### Delete input sample-specific directories if requested
if [[ -v DELETE ]]
then
  printf "\n\n--->Deleting input sample-specific directories...<---\n"
  rm -r $(cat "$SAMPLES")
fi

printf "\n--->JOB COMPLETE<---\n\n"
