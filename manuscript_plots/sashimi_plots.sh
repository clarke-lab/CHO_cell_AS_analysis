#!/usr/bin/env bash
sample1_rep1=../lncRNA_analysis/mapped/REP37_1Aligned.sortedByCoord.out.bam
sample1_rep2=../lncRNA_analysis/mapped/REP37_2Aligned.sortedByCoord.out.bam
sample1_rep3=../lncRNA_analysis/mapped/REP37_3Aligned.sortedByCoord.out.bam
sample1_rep4=../lncRNA_analysis/mapped/REP37_4Aligned.sortedByCoord.out.bam
sample2_rep1=../lncRNA_analysis/mapped/REP31_1Aligned.sortedByCoord.out.bam
sample2_rep2=../lncRNA_analysis/mapped/REP31_2Aligned.sortedByCoord.out.bam
sample2_rep3=../lncRNA_analysis/mapped/REP31_3Aligned.sortedByCoord.out.bam
sample2_rep4=../lncRNA_analysis/mapped/REP31_4Aligned.sortedByCoord.out.bam

output_dir="manuscript_plots/sashimi_plots"
mkdir $output_dir

awk '$1~/^47910$/ {print $0}' rmats/SE.MATS.JC.txt | \
sed 's/chr//g' | \
awk '$3="Dnm1l"' FS="\t" OFS="\t" | \
awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' \
>  $output_dir/Dnm1l.txt

awk '$1~/^8557$/ {print $0}' rmats/SE.MATS.JC.txt | \
sed 's/chr//g' | \
awk '$3="Dnm1l"' FS="\t" OFS="\t" | \
awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' \
>  $output_dir/Gpx8.txt

rmats2sashimiplot \
--b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  \
--b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 \
-e $output_dir/Dnm1l.txt \
-t SE \
--l1 31 \
--l2 37 \
--exon_s 1 \
--intron_s 5 \
-o $output_dir/Dnm1l --group-info data/group.gf --color "#F9B288","#A2D9F9"

rmats2sashimiplot \
--b1 $sample1_rep1,$sample1_rep2,$sample1_rep3,$sample1_rep4  \
--b2 $sample2_rep1,$sample2_rep2,$sample2_rep3,$sample2_rep4 \
-e $output_dir/Gpx8.txt \
-t SE \
--l1 31 \
--l2 37 \
--exon_s 1 \
--intron_s 5 \
-o $output_dir/Gpx8 --group-info data/group.gf --color "#F9B288","#A2D9F9"
