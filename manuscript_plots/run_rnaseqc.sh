#!/usr/bin/env bash
sample1_rep1=../lncRNA_analysis/mapped/REP37_1Aligned.sortedByCoord.out.bam
sample1_rep2=../lncRNA_analysis/mapped/REP37_2Aligned.sortedByCoord.out.bam
sample1_rep3=../lncRNA_analysis/mapped/REP37_3Aligned.sortedByCoord.out.bam
sample1_rep4=../lncRNA_analysis/mapped/REP37_4Aligned.sortedByCoord.out.bam
sample2_rep1=../lncRNA_analysis/mapped/REP31_1Aligned.sortedByCoord.out.bam
sample2_rep2=../lncRNA_analysis/mapped/REP31_2Aligned.sortedByCoord.out.bam
sample2_rep3=../lncRNA_analysis/mapped/REP31_3Aligned.sortedByCoord.out.bam
sample2_rep4=../lncRNA_analysis/mapped/REP31_4Aligned.sortedByCoord.out.bam

awk '$3~/^transcript$/ {print $0}' stringtie/stringtie.gtf |
> manuscript_plots\qc_metrics\stringtie.transcript.gtf


mkdir manuscript_plots/rnaseqc/duplicate_marked

java -jar   /mnt/HDD2/colin/bin/picard.jar MarkDuplicates \
      I="$sample1_rep1" \
      OUTPUT=manuscript_plots/rnaseqc/duplicate_marked/REP37_1_marked_duplicates.bam \
      M=manuscript_plots/rnaseqc/duplicate_marked/REP37_1_marked_dup_metrics.txt

java -jar   /mnt/HDD2/colin/bin/picard.jar MarkDuplicates \
I="$sample1_rep2" O=manuscript_plots/rnaseqc/duplicate_marked/REP37_2_marked_duplicates.bam \
M=manuscript_plots/rnaseqc/duplicate_marked/REP37_2_marked_dup_metrics.txt&

java -jar   /mnt/HDD2/colin/bin/picard.jar MarkDuplicates \
I="$sample1_rep3" \
OUTPUT=manuscript_plots/rnaseqc/duplicate_marked/REP37_3_marked_duplicates.bam \
M=manuscript_plots/rnaseqc/duplicate_marked/REP37_3_marked_dup_metrics.txt&


java -jar   /mnt/HDD2/colin/bin/picard.jar MarkDuplicates \
I="$sample1_rep4" \
OUTPUT=manuscript_plots/rnaseqc/duplicate_marked/REP37_4_marked_duplicates.bam \
M=manuscript_plots/rnaseqc/duplicate_marked/REP37_4_marked_dup_metrics.txt&


java -jar   /mnt/HDD2/colin/bin/picard.jar MarkDuplicates \
I="$sample2_rep1" \
OUTPUT=manuscript_plots/rnaseqc/duplicate_marked/REP31_1_marked_duplicates.bam \
M=manuscript_plots/rnaseqc/duplicate_marked/REP37_1_marked_dup_metrics.txt&


java -jar   /mnt/HDD2/colin/bin/picard.jar MarkDuplicates \
I="$sample2_rep2" \
OUTPUT=manuscript_plots/rnaseqc/duplicate_marked/REP31_2_marked_duplicates.bam \
M=manuscript_plots/rnaseqc/duplicate_marked/REP31_2_marked_dup_metrics.txt&

java -jar   /mnt/HDD2/colin/bin/picard.jar MarkDuplicates \
I="$sample2_rep3" \
OUTPUT=manuscript_plots/rnaseqc/duplicate_marked/REP31_3_marked_duplicates.bam \
M=manuscript_plots/rnaseqc/duplicate_marked/REP31_1_marked_dup_metrics.txt&


java -jar   /mnt/HDD2/colin/bin/picard.jar MarkDuplicates \
I="$sample2_rep4" \
OUTPUT=manuscript_plots/rnaseqc/duplicate_marked/REP31_4_marked_duplicates.bam \
M=manuscript_plots/rnaseqc/duplicate_marked/REP31_1_marked_dup_metrics.txt&


java -jar /mnt/HDD2/colin/bin/picard.jar AddOrReplaceReadGroups \
I=manuscript_plots/rnaseqc/duplicate_marked/REP37_1_marked_duplicates.bam \
O=manuscript_plots/rnaseqc/duplicate_marked/REP37_1_marked_duplicates.rg.bam \
RGSM="REP37_1" \
RGLB=rnaseq \
RGPL=illumina \
RGPU=none \
VALIDATION_STRINGENCY=LENIENT;

java -jar /mnt/HDD2/colin/bin/picard.jar AddOrReplaceReadGroups \
I=manuscript_plots/rnaseqc/duplicate_marked/REP37_2_marked_duplicates.bam \
O=manuscript_plots/rnaseqc/duplicate_marked/REP37_2_marked_duplicates.rg.bam \
RGSM="REP37_2" \
RGLB=rnaseq \
RGPL=illumina \
RGPU=none \
VALIDATION_STRINGENCY=LENIENT;

java -jar /mnt/HDD2/colin/bin/picard.jar AddOrReplaceReadGroups \
I=manuscript_plots/rnaseqc/duplicate_marked/REP37_3_marked_duplicates.bam \
O=manuscript_plots/rnaseqc/duplicate_marked/REP37_3_marked_duplicates.rg.bam \
RGSM="REP37_3" \
RGLB=rnaseq \
RGPL=illumina \
RGPU=none \
VALIDATION_STRINGENCY=LENIENT;

java -jar /mnt/HDD2/colin/bin/picard.jar AddOrReplaceReadGroups \
I=manuscript_plots/rnaseqc/duplicate_marked/REP37_4_marked_duplicates.bam \
O=manuscript_plots/rnaseqc/duplicate_marked/REP37_4_marked_duplicates.rg.bam \
RGSM="REP37_4" \
RGLB=rnaseq \
RGPL=illumina \
RGPU=none \
VALIDATION_STRINGENCY=LENIENT;

java -jar /mnt/HDD2/colin/bin/picard.jar AddOrReplaceReadGroups \
I=manuscript_plots/rnaseqc/duplicate_marked/REP31_1_marked_duplicates.bam \
O=manuscript_plots/rnaseqc/duplicate_marked/REP31_1_marked_duplicates.rg.bam \
RGSM="REP31_1" \
RGLB=rnaseq \
RGPL=illumina \
RGPU=none \
VALIDATION_STRINGENCY=LENIENT;

java -jar /mnt/HDD2/colin/bin/picard.jar AddOrReplaceReadGroups \
I=manuscript_plots/rnaseqc/duplicate_marked/REP31_2_marked_duplicates.bam \
O=manuscript_plots/rnaseqc/duplicate_marked/REP31_2_marked_duplicates.rg.bam \
RGSM="REP31_2" \
RGLB=rnaseq \
RGPL=illumina \
RGPU=none \
VALIDATION_STRINGENCY=LENIENT;

java -jar /mnt/HDD2/colin/bin/picard.jar AddOrReplaceReadGroups \
I=manuscript_plots/rnaseqc/duplicate_marked/REP31_3_marked_duplicates.bam \
O=manuscript_plots/rnaseqc/duplicate_marked/REP31_3_marked_duplicates.rg.bam \
RGSM="REP31_3" \
RGLB=rnaseq \
RGPL=illumina \
RGPU=none \
VALIDATION_STRINGENCY=LENIENT;

java -jar /mnt/HDD2/colin/bin/picard.jar AddOrReplaceReadGroups \
I=manuscript_plots/rnaseqc/duplicate_marked/REP31_4_marked_duplicates.bam \
O=manuscript_plots/rnaseqc/duplicate_marked/REP31_4_marked_duplicates.rg.bam \
RGSM="REP31_4" \
RGLB=rnaseq \
RGPL=illumina \
RGPU=none \
VALIDATION_STRINGENCY=LENIENT;

samtools index manuscript_plots/rnaseqc/duplicate_marked/REP37_1_marked_duplicates.rg.bam&
samtools index manuscript_plots/rnaseqc/duplicate_marked/REP37_2_marked_duplicates.rg.bam&
samtools index manuscript_plots/rnaseqc/duplicate_marked/REP37_3_marked_duplicates.rg.bam&
samtools index manuscript_plots/rnaseqc/duplicate_marked/REP37_4_marked_duplicates.rg.bam&
samtools index manuscript_plots/rnaseqc/duplicate_marked/REP31_1_marked_duplicates.rg.bam&
samtools index manuscript_plots/rnaseqc/duplicate_marked/REP31_2_marked_duplicates.rg.bam&
samtools index manuscript_plots/rnaseqc/duplicate_marked/REP31_3_marked_duplicates.rg.bam&
samtools index manuscript_plots/rnaseqc/duplicate_marked/REP31_4_marked_duplicates.rg.bam&

awk '$3~/^transcript$/ {print $0}' stringtie/stringtie.gtf |
> manuscript_plots/rnaseqc/stringtie.transcript.gtf

awk '$3~/^transcript$/ {print $0}' reference_genome/ensembl_chok1_genome.gtf |
> manuscript_plots/rnaseqc/reference.transcript.gtf


  java -jar /home/colin/bin/RNA-SeQC_v1.1.8.jar \
  -s manuscript_plots/rnaseqc/rnaseq_sample_file.txt \
  -t manuscript_plots/rnaseqc/reference.transcript.gtf \
  -r ../lncRNA_analysis/reference_genome/ensembl_chok1_genome.fa \
  -o manuscript_plots/rnaseqc/reference_annotation


 ~/bin/rnaseqc.v2.3.0.linux \
stringtie/stringtie.gtf \
 $sample1_rep1 \
 --stranded rf \
 output manuscript_plots\qc_metrics\stringtie
