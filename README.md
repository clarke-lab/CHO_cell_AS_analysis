# alt_splicing_analysis

## Prerequisites 
Trimmomatic [modify the code to add to path]
Cutadpat
STAR
R
Samtools

# get the ensembl CHOK1 genome and GTF
bash prepare_genome.sh

#make the STAR index
make_star_index.sh

###data preprocessing
bash preprocess_data.sh total-rna-1 ../data/raw ../data/preprocessed
bash preprocess_data.sh total-rna-2 ../data/raw ../data/preprocessed
bash preprocess_data.sh total-rna-3 ../data/raw ../data/preprocessed
bash preprocess_data.sh total-rna-4 ../data/raw ../data/preprocessed
bash preprocess_data.sh total-rna-5 ../data/raw ../data/preprocessed
bash preprocess_data.sh total-rna-6 ../data/raw ../data/preprocessed
bash preprocess_data.sh total-rna-7 ../data/raw ../data/preprocessed
bash preprocess_data.sh total-rna-8 ../data/raw ../data/preprocessed

###data mapping
'''bash
bash preprocess_data.sh total-rna-1 ../data/preprocessed/paired
bash preprocess_data.sh total-rna-2 ../data/preprocessed/paired
bash preprocess_data.sh total-rna-3 ../data/preprocessed/paired
bash preprocess_data.sh total-rna-4 ../data/preprocessed/paired
bash preprocess_data.sh total-rna-5 ../data/preprocessed/paired
bash preprocess_data.sh total-rna-6 ../data/preprocessed/paired
bash preprocess_data.sh total-rna-7 ../data/preprocessed/paired
'''

###count for DESeq2
bash htseq_count.sh total-rna-1 
bash htseq_count.sh total-rna-2 
bash htseq_count.sh total-rna-3 
bash htseq_count.sh total-rna-4
bash htseq_count.sh total-rna-5
bash htseq_count.sh total-rna-6
bash htseq_count.sh total-rna-7

### string tie assembly
bash stringtie_star.sh total-rna-1
bash stringtie_star.sh total-rna-2 
bash stringtie_star.sh total-rna-3 
bash stringtie_star.sh total-rna-4 
bash stringtie_star.sh total-rna-5
bash stringtie_star.sh total-rna-6 
bash stringtie_star.sh total-rna-7 
bash stringtie_star.sh total-rna-8 

### trim data for rMATs
bash rmats_trim.sh total-rna-1
bash rmats_trim.sh total-rna-2 
bash rmats_trim.sh total-rna-3 
bash rmats_trim.sh total-rna-4 
bash rmats_trim.sh total-rna-5
bash rmats_trim.sh total-rna-6 
bash rmats_trim.sh total-rna-7 
bash rmats_trim.sh total-rna-8








