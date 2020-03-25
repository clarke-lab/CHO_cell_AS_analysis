# CHO cell alternative splicing analysis
## Introduction
These scripts  enable the reproduction of the results of

Tzani et al. Sub physiological temperature induces pervasive alternative splicing in Chinese hamster ovary cells

The experimental design is as follows:

Data availability:

### Dependancies
Trimmomatic 0.36 [modify the code to add to path]
Cutadpat
STAR-2.7.2d
R 3.5
DESeq2_1.22.2
biomaRt_2.38.0
rMats.3.5.2
stringtie
Samtools
## RNASeq data preprocesssing
### Download the data from ENA
This is a simple way to dowload from ENA, for higher speed download use the
Aspera client
Total data download size: ~95G
```bash
mkdir -p data/ena
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/057/SRR10572657/*" -P data/ena || { handle ; error ; }&
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/058/SRR10572658/*" -P data/ena || { handle ; error ; }&
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/059/SRR10572659/*" -P data/ena || { handle ; error ; }&
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/060/SRR10572660/*" -P data/ena || { handle ; error ; }&
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/061/SRR10572661/*" -P data/ena || { handle ; error ; }&
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/062/SRR10572662/*" -P data/ena || { handle ; error ; }&
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/063/SRR10572663/*" -P data/ena || { handle ; error ; }&
wget -q "ftp://ftp.sra.ebi.ac.uk/ena/vol1/fastq/SRR105/064/SRR10572664/*" -P data/ena || { handle ; error ; }
```
### create a list of sample sample_names
```bash
ls ../data/raw/ | sed -n 's/\.fastq.gz$//p' | cut -d_ -f1-2 | uniq > data/sample_names.txt
```
### trim adapter sequences
Adapter trimming for the Illumina TruSeq adapters
```bash
mkdir data/cutadapt
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/trim_adapter.sh -s $sample  -i data/ena -o data/cutadapt&
done
```
### quality trimming
Trimmomatic for removing low quality bases and filtering resulting reads that are  
too short. The script makes two subfolders in the trimmomatic folder for paired and unpaired reads  
```bash
mkdir data/preprocessed
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/trim_quality.sh -s $sample -i data/cutadapt -o data/preprocessed -p 32
done
```
## Read Mapping
### Download the reference genome and NCBI annotation for CHO K1
. The sequence is also prepped for further analysis by retaining only the scaffold ID.
In addition, a complementay annotation file is created from NCBI to help with annotation later
```bash
mkdir reference_genome
./scripts/prepare_genome.sh -v 98 -o reference_genome
```
### make the STAR index
Use the reference genome to build an index for mapping the RNASeq reads
```bash
./scripts/make_star_index.sh -g reference_genome/ensembl_chok1_genome.fa -a reference_genome/ensembl_chok1_genome.gtf -p 32 -d reference_genome&
```

### map to CHOK1 genome
```bash
mkdir data/mapped
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/star_mapping.sh -s $sample -i data/preprocessed/paired -g reference_genome/star_index -o data/mapped -p 32
done
```

## Genome guided assembly
### string tie assembly
```bash
mkdir transcriptome_assembly
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
  ./scripts/run_Stringtie.sh -s $sample -i data/mapped -g reference_genome/ensembl_chok1_genome.gtf -o transcriptome_assembly -p 32
done
```

### merge individual stringtie assemblies and compare to ENSEMBL annotation
```bash
./scripts/stringtie_merge.sh -t transcriptome_assembly -g reference_genome/ensembl_chok1_genome.gtf
```

## Differential gene expression analysis
### count for DESeq2
```bash
mkdir DE_results
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/htseq_count.sh -s $sample -m data/mapped -g reference_genome/ensembl_chok1_genome.gtf -o DE_results&
done
```

### count for DESeq2
```bash
Rscript R/run_deseq2.R "DE_results/counts" "DE_results"
```
## Alternative splicing analysis
### Read trimming
rMats requires that all read are of equal length and each read was trimmed to 125bp
```bash
mkdir data/rmats
cat data/sample_info.txt | cut -f 2 | tail -n 8 | while read sample; do
./scripts/rmats_trim.sh -s $sample -i data/preprocessed/paired  -o data/rmats
done
```

### rMAT analysis
```bash
mkdir AS_results/unfiltered
python2 ~/bin/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py \
--s1 data/ts_replicates.txt \
--s2 data/nts_replicates.txt \
--gtf transcriptome_assembly/stringtie.gtf \
--bi reference_genome/star_index \
--od rMats_output \
-t paired \
--readLength 125 \
--libType fr-firststrand
```

### filter and annotate the rmats results
```bash
Rscript R/filter_rmats_results_v1.R "rMats_output" "AS_results"
```
### Sashimi plot example

```bash
# make a directory to store the plots
output_dir="AS_results/sashimi_plots"
mkdir $output_dir

# make a grouping file
touch $output_dir/grouping.gf 
echo  "NTS: 1-4" >> $output_dir/grouping.gf
echo  "TS: 5-8" >> $output_dir/grouping.gf

# specific the bam files to make the plot
ts_rep1=./data/mapped/SRR10572664Aligned.sortedByCoord.out.bam
ts_rep2=./data/mapped/SRR10572663Aligned.sortedByCoord.out.bam
ts_rep3=./data/mapped/SRR10572662Aligned.sortedByCoord.out.bam
ts_rep4=./data/mapped/SRR10572661Aligned.sortedByCoord.out.bam
nts_rep1=./data/mapped/SRR10572660Aligned.sortedByCoord.out.bam
nts_rep2=./data/mapped/SRR10572659Aligned.sortedByCoord.out.bam
nts_rep3=./data/mapped/SRR10572658Aligned.sortedByCoord.out.bam
nts_rep4=./data/mapped/SRR10572657Aligned.sortedByCoord.out.bam

# select the rMats event to plot from the UNFILTERED out
awk '$1~/^14775$/ {print $0}' rMats_output/SE.MATS.JC.txt | \
sed 's/chr//g' | \
awk '$3="Dnm1l"' FS="\t" OFS="\t" | \
awk -F $'\t' ' { t = $21; $21 = $22; $22 = t; print; } ' OFS=$'\t' \
>  $output_dir/Dnm1l.txt

# generate the sashimi plot
rmats2sashimiplot \
--b2 $ts_rep1,$ts_rep2,$ts_rep3,$ts_rep4  \
--b1 $nts_rep1,$nts_rep2,$nts_rep3,$nts_rep4 \
-e $output_dir/Dnm1l.txt \
-t SE \
--l1 31 \
--l2 37 \
--exon_s 1 \
--intron_s 5 \
-o $output_dir/Dnm1l --group-info $output_dir/grouping.gf --color "#F9B288","#A2D9F9"
```