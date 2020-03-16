#!/usr/bin/env bash
#### get the chinese hamster reference genome from ENSEMBL
#### get additional annotations from NCBI
#### inputs are: 1) output directory and 2) ensembl version
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

while getopts v:o: option
  do
    case "${option}"
      in
      v) ENSEMBL_VER=${OPTARG};;
      o) GENOME_DIR=${OPTARG};;
    esac
done

mkdir -p $GENOME_DIR

# get the ENSEMBL CHOK1 reference genome and GTF file
wget ftp://ftp.ensembl.org/pub/release-$ENSEMBL_VER/fasta/cricetulus_griseus_crigri/dna/Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa.gz \
-P $GENOME_DIR
wget ftp://ftp.ensembl.org/pub/release-$ENSEMBL_VER/gtf/cricetulus_griseus_crigri/Cricetulus_griseus_crigri.CriGri_1.0.98.gtf.gz \
-P $GENOME_DIR

gunzip $GENOME_DIR/*.gz

# retain only scaffold ID in the fasta header
raw_reference="Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa"
name="ensembl_chok1_genome"
sed '/^>/ s/ .*//' $GENOME_DIR/$raw_reference > $GENOME_DIR/$name.fa
mv $GENOME_DIR/Cricetulus_griseus_crigri.CriGri_1.0.98.gtf $GENOME_DIR/$name.gtf
rm $GENOME_DIR/$raw_reference

# get NCBI CHOK1 Entrez annotations and create an annotation file
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz \
-P $GENOME_DIR

gunzip $GENOME_DIR/*.gz
grep '^\<10029\>' $GENOME_DIR/gene_info > $GENOME_DIR/chok1_ncbi_ids.txt
rm $GENOME_DIR/gene_info
#list on ensembl protein coding genes
#grep 'protein_coding' $GENOMEDIR/$name.gtf | awk {'print $10'} | uniq | tr -d \" | tr -d \; > $GENOMEDIR/protein.coding.genes.list


#list of ensembl other ncRNAs
#egrep 'miRNA|Mt_rRNA|Mt_tRNA|processed_pseudogene|pseudogene|ribozyme|rRNA|scaRNA|snoRNA|snRNA|sRNA' $GENOMEDIR/$name.gtf \
#| awk {'print $10'} | uniq | tr -d \" | tr -d \; > $TRANSCRIPTDIR/other.noncoding.genes.list

# make lists for lncrNAs annotated in ensembl
# grep 'lincRNA' $GENOMEDIR/$name.gtf | awk {'print $10'} | uniq | tr -d \" | tr -d \; > $GENOMEDIR/ensembl_lncrna_gene.list
# grep 'lincRNA' $GENOMEDIR/$name.gtf | awk {'print $14'} | grep -v 'ensembl' | uniq | tr -d \" | tr -d \; > $GENOMEDIR/ensembl_lncrna_transcript.list
