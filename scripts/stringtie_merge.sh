#!/usr/bin/env bash

if (($# == 0)); then
        echo "Usage:"
        echo "-t = assembled transcript directory"
        echo "-g = path to reference annotation"
        exit 2
fi
while getopts t:g: option
  do
    case "${option}"
      in
      t) TRANSCRIPTDIR=${OPTARG};;
      g) GTF=${OPTARG};;
    esac
done


readlink -f $TRANSCRIPTDIR/individual_gtfs/*.gtf >> $TRANSCRIPTDIR/mergelist.txt

/mnt/HDD2/colin/bin/stringtie-2.0.3.Linux_x86_64/stringtie \
--merge $TRANSCRIPTDIR/mergelist.txt \
-o $TRANSCRIPTDIR/stringtie_original.gtf \
-G $GTF \
-f 0.1  \
-c 10

# create a file liniking stringtie ID tO ENSEMBL geNE ID
grep -wFf reference_genome/protein.coding.genes.list $TRANSCRIPTDIR/stringtie_original.gtf | \
grep -v exon | awk '{print $10, $NF}' | uniq | tr -d \" | tr -d \; > $TRANSCRIPTDIR/stringtie_ensembl_gene_mapping.list

# append ensembl gene ids to MSTRG GTF
perl scripts/mstrg_prep.pl $TRANSCRIPTDIR/stringtie_original.gtf > $TRANSCRIPTDIR/stringtie_merged.appended.gtf

# find instances where stringtie has asemebled transcripts from 2 or more overlaping loci and created a new "gene".
# The final field of the GTF file will contain an MSTRG ID not an ENS ID
grep 'MSTRG.*|ENSCGRG.*|ENSC.*' $TRANSCRIPTDIR/stringtie_merged.appended.gtf | \
grep '\<transcript\>' | awk '$NF ~/MSTRG/ {print $NF}'  > $TRANSCRIPTDIR/removed.overlapped.MSTRG.transcripts

# remove assembled transcripts spanning two or more sense overlapping genes transcripts
grep -v -F -f $TRANSCRIPTDIR/removed.overlapped.MSTRG.transcripts $TRANSCRIPTDIR/stringtie_merged.appended.gtf > $TRANSCRIPTDIR/stringtie_merged.appended.fp.filtered.gtf

# remove transcripts without strand
awk '$7 != "." {print}' $TRANSCRIPTDIR/stringtie_merged.appended.fp.filtered.gtf > $TRANSCRIPTDIR/stringtie.gtf

gffcompare \
-o $TRANSCRIPTDIR/gffcompare \
-r $GTF $TRANSCRIPTDIR/stringtie.gtf
