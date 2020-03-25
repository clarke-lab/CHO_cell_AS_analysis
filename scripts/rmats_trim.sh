#!/usr/bin/env bash
#### Clip the RNAseq reads to a standard read length
#### inputs are: 1) sample ID and 2) untrimmed FASTQ directory 3) output directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2020

if (($# == 0)); then
        echo "Usage:"
        echo "-s = sample ID"
        echo "-i = input directory"
        echo "-o = path for output BAMs"
        exit 2
fi
while getopts s:i:o: option
  do
    case "${option}"
      in
      s) SAMPLE_ID=${OPTARG};;
      i) IN_DIR=${OPTARG};;
      o) OUT_DIR=${OPTARG};;
    esac
done

mkdir -p $OUT_DIR/trimmed
mkdir -p $OUT_DIR/untrimmed

# FORWARD READ
cp $IN_DIR/"$SAMPLE_ID"_1.fastq.gz $OUT_DIR/untrimmed

gunzip $OUT_DIR/untrimmed/"$SAMPLE_ID"_1.fastq.gz

python ~/bin/rMATS.3.2.5/bin/trimFastq.py  \
$OUT_DIR/untrimmed/"$SAMPLE_ID"_1.fastq \
$OUT_DIR/trimmed/"$SAMPLE_ID"_1.fastq 125
rm $OUT_DIR/untrimmed/"$SAMPLE_ID"_1.fastq

# REVERSE READ
cp $IN_DIR/"$SAMPLE_ID"_2.fastq.gz $OUT_DIR/untrimmed

gunzip $OUT_DIR/untrimmed/"$SAMPLE_ID"_2.fastq.gz

python ~/bin/rMATS.3.2.5/bin/trimFastq.py  \
$OUT_DIR/untrimmed/"$SAMPLE_ID"_2.fastq \
$OUT_DIR/trimmed/"$SAMPLE_ID"_2.fastq 125
rm $OUT_DIR/untrimmed/"$SAMPLE_ID"_2.fastq
