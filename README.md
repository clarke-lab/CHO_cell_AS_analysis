# alt_splicing_analysis

These scripts replicate the results of the following manuscript

## Installation

### Dependancies




## Prerequisites
Trimmomatic [modify the code to add to path]
Cutadpat
STAR
R
Samtools

### get the ensembl CHOK1 genome and GTF
```bash
./scripts/prepare_genome.sh -v 98 -o reference_genome
```
#make the STAR index

```bash
./scripts/make_star_index.sh -g reference_genome/ensembl_chok1_genome.fa -a reference_genome/ensembl_chok1_genome.gtf -p 32 -d reference_genome
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

### map to CHOK1 genome
```bash
cat data/sample_names.txt | while read sample; do
./scripts/star_mapping.sh -s $sa mple -i ../data/preprocessed/paired -g reference_genome/star_index -o mapped -p 32
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
