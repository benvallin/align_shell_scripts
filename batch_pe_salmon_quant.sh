#!/bin/bash

### Define help message
HELP='
NAME
  batch_pe_salmon_quant.sh - run salmon quant iteratively on a set of PE samples 

SYNOPSIS
  batch_pe_salmon_quant.sh -s -i -l -1 -2 -o [-p -g -S -d -D -w -b]
  batch_pe_salmon_quant.sh -h

DESCRIPTION
  Run salmon quant on a set of PE samples and optionally produce coordinate-sorted BAM files of selective-alignment results.
  -s  File listing the samples to operate on.
      Each line should contain the path to a sample-specific directory containing pair-end reads fastq files.
      Path can be absolute or relative to the current directory.
  -i  Path to Salmon index (passed as -i to salmon quant). 
      Path can be absolute or relative to the current directory.
  -l  Format string describing the library type (passed as -l to salmon quant).
  -1  Regexp matching the name of fastq files containing the #1 mates (passed as -1 to salmon quant).
  -2  Regexp matching the name of fastq files containing the #2 mates (passed as -2 to salmon quant).
  -o  Path to output quantification directory.
      Results are written to sample-specific subdirectories within the provided directory.
      Path can be absolute or relative to the current directory.
  -p  The number of threads to use concurrently (passed as -p to salmon quant). Default is 8.
  -g  Perform fragment GC bias correction (passed as --gcBias to salmon quant).
  -S  Perform sequence-specific bias correction (passed as --seqBias to salmon quant).
  -d  Allow dovetailing mappings (passed as --allowDovetail to salmon quant).
  -D  Dump the simple equivalence class counts that were computed during mapping or alignment (passed as --dumpEq to salmon quant).
  -w  Write the names of un-mapped reads to the file unmapped_names.txt in the auxiliary directory (passed as --writeUnmappedNames to salmon quant).
  -b  Produce coordinate-sorted BAM files of selective-alignment results (pass --writeUnmappedNames to salmon quant and pipe to samtools view/sort).
  -h  Print this help.\n
'

### Parse option arguments
while getopts ":hs:i:l:1:2:o:p:gSdDwb" opt
do
  case $opt in
    h) printf "$HELP"
    exit 0
    ;;
    s) SAMPLES="$OPTARG"
    ;;
    i) INDEX="$OPTARG"
    ;;
    l) LIBTYPE="$OPTARG"
    ;;
    1) MATES1="$OPTARG"
    ;;
    2) MATES2="$OPTARG"
    ;;
    o) OUTDIR="${OPTARG%/}"
    ;;
    p) THREADS="$OPTARG"
    ;;
    g) GCBIAS="--gcBias"
    ;;
    S) SEQBIAS="--seqBias"
    ;;
    d) ALLOWDOVETAIL="--allowDovetail"
    ;;
    D) DUMPEQ="--dumpEq"
    ;;
    w) WRITEUNMAPPEDNAMES="--writeUnmappedNames"
    ;;
    b) BAM="--writeMappings"
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

# => Set n threads to default value if not provided
if ! [[ -v THREADS ]]
then 
  THREADS=8
fi

### Construct Salmon's arguments
SALMON_ARGS=()

# => Add systematic arguments to Salmon's arguments
SALMON_ARGS+=(--index $INDEX --libType $LIBTYPE --validateMappings --threads $THREADS)

# => Add --gcBias, --seqBias, --allowDovetail, --dumpEq and --writeUnmappedNames to Salmon's arguments if provided 
OPTIONAL_ARGS="$GCBIAS $SEQBIAS $ALLOWDOVETAIL $DUMPEQ $WRITEUNMAPPEDNAMES"

for i in "$OPTIONAL_ARGS"
do SALMON_ARGS+=($i)
done

### Record samples details
N_SAMPLES=$(cat "$SAMPLES" | wc -l)
CURRENT_SAMPLE=0

printf "\n\n--->$N_SAMPLES SAMPLES TO PROCESS<---"

### Run salmon quant
for i in $(cat "$SAMPLES")
do 
  CURRENT_SAMPLE=$((CURRENT_SAMPLE+1))
  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES: PROCESSING SAMPLE $(basename $i)<---\n"

  # => If user requested coordinate-sorted BAM files 
  if [[ -v BAM ]]
  then
    printf "\n--->Selective-alignment results will be written to coordinate-sorted BAM file<---\n\n" 
    salmon quant "${SALMON_ARGS[@]}" \
    --mates1 <(gunzip -c $(find "$i" -name "$MATES1" | sort)) \
    --mates2 <(gunzip -c $(find "$i" -name "$MATES2" | sort)) \
    --output "$OUTDIR"/$(basename "$i")/ \
    --writeMappings | \
    samtools view -Shu -@ "$THREADS" - | \
    samtools sort -@ "$THREADS" -o "$OUTDIR"/$(basename "$i")/$(basename "$i")_writeMappings_output_coordinates_sorted.bam -

  # => If user did not request coordinate-sorted BAM files 
  else
    printf "\n"
    salmon quant "${SALMON_ARGS[@]}" \
    --mates1 <(gunzip -c $(find "$i" -name "$MATES1" | sort)) \
    --mates2 <(gunzip -c $(find "$i" -name "$MATES2" | sort)) \
    --output "$OUTDIR"/$(basename "$i")/
  fi

  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES COMPLETE<---\n"
done

printf "\n\n--->ALL JOBS COMPLETE<---\n\n"
