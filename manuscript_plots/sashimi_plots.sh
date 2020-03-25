#!/usr/bin/env bash
sample1_rep1=./data/mapped/SRR10572664Aligned.sortedByCoord.out.bam
sample1_rep2=./data/mapped/SRR10572663Aligned.sortedByCoord.out.bam
sample1_rep3=./data/mapped/SRR10572662Aligned.sortedByCoord.out.bam
sample1_rep4=./data/mapped/SRR10572661Aligned.sortedByCoord.out.bam
sample2_rep1=./data/mapped/SRR10572660Aligned.sortedByCoord.out.bam
sample2_rep2=./data/mapped/SRR10572659Aligned.sortedByCoord.out.bam
sample2_rep3=./data/mapped/SRR10572658Aligned.sortedByCoord.out.bam
sample2_rep4=./data/mapped/SRR10572657Aligned.sortedByCoord.out.bam

output_dir="AS_results/sashimi_plots"
mkdir $output_dir

touch $output_dir/grouping.gf 
echo  "NTS: 1-4" >> $output_dir/grouping.gf
echo  "TS: 5-8" >> $output_dir/grouping.gf


awk '$1~/^14775$/ {print $0}' rMats_output/SE.MATS.JC.txt | \
sed 's/chr//g' | \
awk '$3="Dnm1l"' FS="\t" OFS="\t" | \
awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' \
>  $output_dir/Dnm1l.txt

awk '$1~/^51352$/ {print $0}' rMats_output/SE.MATS.JC.txt | \
sed 's/chr//g' | \
awk '$3="Gpx8"' FS="\t" OFS="\t" | \
awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' \
>  $output_dir/Gpx8.txt

rmats2sashimiplot \
--b2 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  \
--b1 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 \
-e $output_dir/Dnm1l.txt \
-t SE \
--l1 31 \
--l2 37 \
--exon_s 1 \
--intron_s 5 \
-o $output_dir/Dnm1l --group-info $output_dir/grouping.gf --color "#F9B288","#A2D9F9"

rmats2sashimiplot \
--b2 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  \
--b1 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 \
-e $output_dir/Gpx8.txt \
-t SE \
--l1 31 \
--l2 37 \
--exon_s 1 \
--intron_s 5 \
-o $output_dir/Gpx8 --group-info $output_dir/grouping.gf --color "#F9B288","#A2D9F9"

