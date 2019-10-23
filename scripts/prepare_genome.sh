#!/usr/bin/env bash
while getopts v:o: option
  do
    case "${option}"
      in
      v) ENSEMBLV=${OPTARG};;
      o) GENOMEDIR=${OPTARG};;
    esac
done

mkdir -p $GENOMEDIR

# get the ENSEMBL CHOK1 reference genome and GTF file
wget ftp://ftp.ensembl.org/pub/release-$ENSEMBLV/fasta/cricetulus_griseus_crigri/dna/Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa.gz \
-P $GENOMEDIR
wget ftp://ftp.ensembl.org/pub/release-$ENSEMBLV/gtf/cricetulus_griseus_crigri/Cricetulus_griseus_crigri.CriGri_1.0.98.gtf.gz \
-P $GENOMEDIR

gunzip $GENOMEDIR/*.gz

# retain only scaffold ID in the fasta header
raw_reference="Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa"
name="ensembl_chok1_genome"
sed '/^>/ s/ .*//' $GENOMEDIR/$raw_reference > $GENOMEDIR/$name.fa
mv $GENOMEDIR/Cricetulus_griseus_crigri.CriGri_1.0.98.gtf $GENOMEDIR/$name.gtf
rm $GENOMEDIR/$raw_reference

# get NCBI CHOK1 Entrez annotations
wget ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz \
-P $GENOMEDIR

gunzip $GENOMEDIR/*.gz
grep '^\<10029\>' $GENOMEDIR/gene_info > $GENOMEDIR/chok1_ncbi_ids.txt
rm $GENOMEDIR/gene_info


#list on ensembl protein coding genes
grep 'protein_coding' $GENOMEDIR/$name.gtf | awk {'print $10'} | uniq | tr -d \" | tr -d \; > $GENOMEDIR/protein.coding.genes.list


#list of ensembl other ncRNAs
egrep 'miRNA|Mt_rRNA|Mt_tRNA|processed_pseudogene|pseudogene|ribozyme|rRNA|scaRNA|snoRNA|snRNA|sRNA' $GENOMEDIR/$name.gtf \
| awk {'print $10'} | uniq | tr -d \" | tr -d \; > $TRANSCRIPTDIR/other.noncoding.genes.list

# make lists for lncrNAs annotated in ensembl
grep 'lincRNA' $GENOMEDIR/$name.gtf | awk {'print $10'} | uniq | tr -d \" | tr -d \; > $GENOMEDIR/ensembl_lncrna_gene.list
grep 'lincRNA' $GENOMEDIR/$name.gtf | awk {'print $14'} | grep -v 'ensembl' | uniq | tr -d \" | tr -d \; > $GENOMEDIR/ensembl_lncrna_transcript.list
