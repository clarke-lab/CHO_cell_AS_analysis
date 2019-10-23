#!/bin/sh
raw_dir=../data/raw
out_dir=manuscript_plots/read_counts
cat ../data/sample_names.txt | while read sample; do
num_forward_reads=$(zcat $raw_dir/"$sample"_R1.fastq.gz | echo $((`wc -l`/4)))
num_reverse_reads=$(zcat $raw_dir/"$sample"_R2.fastq.gz | echo $((`wc -l`/4)))
echo $sample'\t'$num_forward_reads'\t'$num_reverse_reads >> $out_dir/$sample.raw.data.counts
done

preprocessed_dir=../data/preprocessed
out_dir=manuscript_plots/read_counts
cat ../data/sample_names.txt | while read sample; do
num_forward_reads=$(zcat $preprocessed_dir/"$sample"_R1.fastq.gz | echo $((`wc -l`/4)))
num_reverse_reads=$(zcat $preprocessed_dir/"$sample"_R2.fastq.gz | echo $((`wc -l`/4)))
echo $sample'\t'$num_forward_reads'\t'$num_reverse_reads >> $out_dir/$sample.preprocessed_dir.data.counts
done
  
