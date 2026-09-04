#!/bin/bash

### Define help message
HELP='
NAME
  batch_pe_gatk_star_align.sh - run STAR iteratively on a set of PE samples for variant calling (GATK pipeline)

SYNOPSIS
  batch_pe_gatk_star_align.sh -s -1 -2 -g -o [-j -d -p -b]
  batch_pe_gatk_star_align.sh -h

DESCRIPTION
  Run STAR --runMode alignReads on a set of PE samples for variant calling purpose and as instructed in GATK best practices.
  The following read group information is inserted in the output bam file: ID:flowcellID_laneID_sampleID, SM:sampleID, LB:sampleID, PL:ILLUMINA.
  -s  File listing the samples to operate on.
      Each line should contain the path to a sample-specific directory containing pair-end reads fastq files.
      Path can be absolute or relative to the current directory.
  -1  find regexp matching the name of fastq files containing the #1 mates.
  -2  find regexp matching the name of fastq files containing the #2 mates.
  -g  Path to the directory where genome files are stored (passed as --genomeDir to STAR). 
      Path can be absolute or relative to the current directory.
  -o  Path to output quantification directory (used to construct STAR --outFileNamePrefix).
      Results are written to sample-specific subdirectories within the provided directory.
      Path can be absolute or relative to the current directory.
  -j  Length of the donor/acceptor sequence on each side of the junctions, ideally = (mate length - 1) (passed as --sjdbOverhang to STAR). Default is 100.
  -d  Command that generates FASTA/FASTQ text from each input file and send it to stdout (passed as --readFilesCommand to STAR). Default is \"pigz -dc\".
  -p  The number of threads to use concurrently (passed as --runThreadN to STAR). Default is 8.
  -b  The number of genome bins for coordinate-sorting. Default is 50.
  -h  Print this help.\n
'

### Parse option arguments
while getopts ":hs:1:2:g:o:j:d:p:b:" opt
do
  case $opt in
    h) printf "$HELP"
    exit 0
    ;;
    s) SAMPLES="$OPTARG"
    ;;
    1) MATES1="$OPTARG"
    ;;
    2) MATES2="$OPTARG"
    ;;
    g) GENOMEDIR="$OPTARG"
    ;;
    o) OUTDIR="${OPTARG%/}"
    ;;
    j) SJDBOVERHANG="$OPTARG"
    ;;
    d) RFCOMMAND="$OPTARG"
    ;;
    p) THREADS="$OPTARG"
    ;;
    b) GENOMEBINS="$OPTARG"
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
 
# => Set length of donor/acceptor sequence on each side of the junctions to default value if not provided
if ! [[ -v SJDBOVERHANG ]]
then 
  SJDBOVERHANG=100
fi

# => Set read file command to default value if not provided
if ! [[ -v RFCOMMAND ]]
then 
  RFCOMMAND="pigz -dc"
fi

# => Set n threads to default value if not provided
if ! [[ -v THREADS ]]
then 
  THREADS=8
fi

# => Set number of genome bins for coordinate-sorting to default value if not provided
if ! [[ -v GENOMEBINS ]]
then 
  GENOMEBINS=50
fi

### Perform alignment
N_SAMPLES=$(cat "$SAMPLES" | wc -l)
CURRENT_SAMPLE=0

printf "\n\n--->$N_SAMPLES SAMPLES TO PROCESS<---"

for i in $(cat "$SAMPLES")
do 
  CURRENT_SAMPLE=$((CURRENT_SAMPLE+1))
  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES: PROCESSING SAMPLE $(basename $i)<---\n"
  
  # => Define input fastq files
  R1FILES=$(find "$i" -name "$MATES1" | sort | sed ':a;N;$!ba;s/\n/,/g') 
  R2FILES=$(find "$i" -name "$MATES2" | sort | sed ':a;N;$!ba;s/\n/,/g') 
  
  # => Construct RG lines to be insterted in output bam file
  RGSM="$(basename $i)"
  RGLINES=()
  for j in $(find "$i" -name "$MATES1" | sort)
  do 
    SEQINFO=$($RFCOMMAND "$j" | head -1 | cut -d " " -f 1)
    FLOWCELL=$(echo "$SEQINFO" | cut -d ":" -f 3)
    LANE=$(echo "$SEQINFO" | cut -d ":" -f 4)
    RGID="$FLOWCELL"_"$LANE"_"$RGSM"
    RGLINE="ID:$RGID  SM:$RGSM  LB:$RGSM  PL:ILLUMINA"
    if [ $(echo $RGLINES | wc -w) -eq 0 ]
    then
      RGLINES+=($RGLINE)
    else
      RGLINES+=(" , $RGLINE")
    fi
  done
  RGLINES="${RGLINES[@]}"
  printf "\nR1FILES=$R1FILES\nR2FILES=$R2FILES\n\nRGLINES=$RGLINES\n\n"
  
  # => Run STAR 
  STAR --runMode alignReads \
  --genomeDir "$GENOMEDIR" \
  --runThreadN "$THREADS" \
  --readFilesIn "$R1FILES" "$R2FILES" \
  --readFilesCommand "$RFCOMMAND" \
  --sjdbOverhang "$SJDBOVERHANG" \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic \
  --limitOutSJcollapsed 1000000 \
  --outFileNamePrefix "$OUTDIR"/"$(basename $i)"/"$(basename $i)"_ \
  --outSAMattrRGline $RGLINES \
  --outBAMsortingBinsN "$GENOMEBINS"
  
  printf "\n\n--->JOB $CURRENT_SAMPLE/$N_SAMPLES COMPLETE<---\n"
done

printf "\n\n--->ALL JOBS COMPLETE<---\n\n"
