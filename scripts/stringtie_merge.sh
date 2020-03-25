#!/usr/bin/env bash
#### Merge the individual transcript assemblies
#### inputs are: 1) assembled transcript parent directory and 2) the refernece GTF
#### 3) reference genome directory
#### Written by NIBRT: colin.clarke@nibrt.ie 12-2019

if (($# == 0)); then
        echo "Usage:"
        echo "-t = assembled transcript directory"
        echo "-g = path to reference annotation"
        echo "-r reference_genome_dir"
        exit 2
fi
while getopts t:g:r: option
  do
    case "${option}"
      in
      t) TRANSCRIPT_DIR=${OPTARG};;
      g) GTF=${OPTARG};;
      r) REF_DIR=${OPTARG};;
    esac
done

readlink -f $TRANSCRIPT_DIR/individual_gtfs/*.gtf >> $TRANSCRIPT_DIR/mergelist.txt

stringtie \
--merge $TRANSCRIPT_DIR/mergelist.txt \
-o $TRANSCRIPT_DIR/stringtie_original.gtf \
-G $GTF \
-f 0.1  \
-c 10

# create a file liniking stringtie ID tO ENSEMBL geNE ID
grep -wFf $REF_DIR/protein.coding.genes.list $TRANSCRIPT_DIR/stringtie_original.gtf | \
grep -v exon | awk '{print $10, $NF}' | uniq | tr -d \" | tr -d \; > $TRANSCRIPT_DIR/stringtie_ensembl_gene_mapping.list

# append ensembl gene ids to MSTRG GTF
perl scripts/mstrg_prep.pl $TRANSCRIPT_DIR/stringtie_original.gtf > $TRANSCRIPT_DIR/stringtie_merged.appended.gtf

# find instances where stringtie has asemebled transcripts from 2 or more overlaping loci and created a new "gene".
# The final field of the GTF file will contain an MSTRG ID not an ENS ID
grep 'MSTRG.*|ENSCGRG.*|ENSC.*' $TRANSCRIPT_DIR/stringtie_merged.appended.gtf | \
grep '\<transcript\>' | awk '$NF ~/MSTRG/ {print $NF}'  > $TRANSCRIPT_DIR/removed.overlapped.MSTRG.transcripts

# remove assembled transcripts spanning two or more sense overlapping genes transcripts
grep -v -F -f $TRANSCRIPT_DIR/removed.overlapped.MSTRG.transcripts $TRANSCRIPT_DIR/stringtie_merged.appended.gtf > $TRANSCRIPT_DIR/stringtie_merged.appended.fp.filtered.gtf

# remove transcripts without strand
awk '$7 != "." {print}' $TRANSCRIPT_DIR/stringtie_merged.appended.fp.filtered.gtf > $TRANSCRIPT_DIR/stringtie.gtf

gffcompare \
-o $TRANSCRIPT_DIR/gffcompare \
-r $GTF $TRANSCRIPT_DIR/stringtie.gtf

# END
