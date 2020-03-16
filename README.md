# alt_splicing_analysis

## Introduction
These scripts  enable the reproduction of the results presented in the following manuscript

Tzani et al. Sub physiological temperature induces pervasive alternative splicing in Chinese hamster ovary cells

The experimental design is as follows:

Data availability:

### Dependancies
Trimmomatic [modify the code to add to path]
Cutadpat
STAR
R (>=3.5)
rMats.3.5.2
stringtie 
Samtools

### Download the data from ENA
This is a simple way to dowload from ENA, for higher speed download use the
Aspera client
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
```bash
cat data/sample_names.txt | while read sample; do
./scripts/trim_adapter.sh -s $sample  -i ../data/raw -o ../data/preprocessed/cutadapt&
done
```

### quality trimming
```bash
cat data/sample_names.txt | while read sample; do
./scripts/trim_quality.sh -s $sample -i ../data/preprocessed/cutadapt -o../data/preprocessed
done
```

## Mapping the RNAseq data
### Download the reference genome and NCBI annotation for CHO K1
. The sequence is also prepped for further analysis by retaining only the scaffold ID.
In addition, a complementay annotation file is created from NCBI to help with annotation later
```bash
./scripts/prepare_genome.sh -v 98 -o reference_genome
```
### make the STAR index
Use the reference genome to build an index for mapping the RNASeq reads
```bash
./scripts/make_star_index.sh -g reference_genome/ensembl_chok1_genome.fa -a reference_genome/ensembl_chok1_genome.gtf -p 32 -d reference_genome
```

### map to CHOK1 genome
```bash
cat data/sample_names.txt | while read sample; do
./scripts/star_mapping.sh -s $sample -i ../data/preprocessed/paired -g reference_genome/star_index -o mapped -p 32
done
```

### count for DESeq2
```bash
mkdir differential_expression
cat ../data/sample_names.txt | while read sample; do
./scripts/htseq_count.sh -s $sample -m ../lncRNA_analysis/mapped -g reference_genome/ensembl_chok1_genome.gtf -o differential_expression/counts&
done
```

### count for DESeq2
```bash
Rscript R/run_DEseq2.R \
  "differential_expression/counts" \
  "differential_expression/DESeq2_results"
```

### string tie assembly
```bash
cat data/sample_names.txt | while read sample; do
./scripts/run_Stringtie.sh -s $sample -i mapped -g reference_genome/ensembl_chok1_genome.gtf -o stringtie -p 32
done
```

### merge individual stringtie assemblies and compare to ENSEMBL annotation
```bash
./scripts/stringtie_merge.sh -t stringtie -g reference_genome/ensembl_chok1_genome.gtf
```



### trim data to fixed 125bp length for rMATs
```bash
cat data/sample_names.txt | while read sample; do
bash rmats_trim.sh REP37_1  ../data/preprocessed/paired
bash rmats_trim.sh REP37_2  ../data/preprocessed/paired
bash rmats_trim.sh REP37_3  ../data/preprocessed/paired
bash rmats_trim.sh REP37_4  ../data/preprocessed/paired
bash rmats_trim.sh REP31_1  ../data/preprocessed/paired
bash rmats_trim.sh REP31_2  ../data/preprocessed/paired
bash rmats_trim.sh REP31_3  ../data/preprocessed/paired
bash rmats_trim.sh REP31_4  ../data/preprocessed/paired
done
```


### rMAT analysis
```bash
python ~/bin/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py \
--s1 data/sample_reps31.txt \
--s2 data/sample_reps37.txt \
--gtf stringtie/stringtie.gtf \
--bi ../lncRNA_analysis/reference_genome/star_index \
--od rmats \
-t paired \
--readLength 125 \
--libType fr-firststrand
```



### filter and annotate the rmats results
```bash
Rscript rscripts/filter_rmats_results.R
```
