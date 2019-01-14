~/bin/TransDecoder-TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl stringtie_output/rmats_stringtie.gtf reference_genome/ensembl_chok1_genome.fa > transdecoder_analysis/stringtie_transcripts.fasta
~/bin/TransDecoder-TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl stringtie_output/rmats_stringtie.gtf > transdecoder_analysis/transcripts.gff3

TransDecoder.LongOrfs -t transdecoder_analysis/stringtie_transcripts.fasta

mkdir pfam_db
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 
gunzip Pfam-A.hmm.gz 

hmmpress pfam_db/Pfam-A.hmm
hmmscan --cpu 28 --domtblout transdecoder_analysis/pfam.domtblout pfam_db/Pfam-A.hmm transdecoder_analysis/stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep -f transdecoder_analysis/pfam.std.out


mkdir transdecoder_analysis/uniref90
wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz 
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz


makeblastdb -in uniprot_sprot.fasta -dbtype prot

blastp -query transdecoder_analysis/stringtie_transcripts.fasta.transdecoder_dir/longest_orfs.pep -db transdecoder_analysis/swissprot/uniprot_sprot.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 32 > transdecoder_analysis/blastp.outfmt6
    
TransDecoder.Predict -t stringtie_transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --single_best_only 



grep -E 'ENSCGRT00000020673|MSTRG.4741.1|MSTRG.4741.2|ENSCGRT00000020672' transdecoder_analysis/stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > qpcr_isoform_proteins/rnf_transcript_ids.txt
seqtk subseq transdecoder_analysis/stringtie_transcripts.fasta.transdecoder.pep qpcr_isoform_proteins/rnf_transcript_ids.txt > qpcr_isoform_proteins/rnf_isoform.proteins.fasta 

~/git/genomeGTFtools/pfampipeline.py rnf_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001 --single_best_only


grep -E 'ENSCGRT00000022494|ENSCGRT00000022493|MSTRG.23141.3|MSTRG.23141.4|MSTRG.23141.5|ENSCGRT00000022496|ENSCGRT00000022495' stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > dnm1l_transcript_ids.txt
seqtk subseq stringtie_transcripts.fasta.transdecoder.pep dnm1l_transcript_ids.txt > dnm1l_isoforms.fasta

~/git/genomeGTFtools/pfampipeline.py dnm1l_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001

grep -E 'MSTRG.21801.2|ENSCGRT00000004998' stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > gpx8_transcript_ids.txt
seqtk subseq stringtie_transcripts.fasta.transdecoder.pep gpx8_transcript_ids.txt > gpx8_isoforms.fasta

~/git/genomeGTFtools/pfampipeline.py gpx8_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001
 

MSTRG.24858.1|MSTRG.24858.2|ENSCGRT00000023840|ENSCGRT00000023837|ENSCGRT00000023838|ENSCGRT00000023839|MSTRG.24858.7|MSTRG.24858.8|ENSCGRT00000023841|ENSCGRT00000023842

grep -E 'MSTRG.24858.1|MSTRG.24858.2|ENSCGRT00000023840|ENSCGRT00000023837|ENSCGRT00000023838|ENSCGRT00000023839|MSTRG.24858.7|MSTRG.24858.8|ENSCGRT00000023841|ENSCGRT00000023842' stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > mff_transcript_ids.txt
seqtk subseq stringtie_transcripts.fasta.transdecoder.pep mff_transcript_ids.txt > mff_isoforms.fasta

~/git/genomeGTFtools/pfampipeline.py mff_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001


grep -E 'ENSCGRT00000015631|MSTRG.14374.2|MSTRG.14374.3|ENSCGRT00000015632' stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > ubr2_transcript_ids.txt
seqtk subseq stringtie_transcripts.fasta.transdecoder.pep ubr2_transcript_ids.txt > ubr2_isoforms.fasta

~/git/genomeGTFtools/pfampipeline.py ubr2_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001


grep -E 'ENSCGRT00000024393|MSTRG.6415.1' stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > slirp_transcript_ids.txt
seqtk subseq stringtie_transcripts.fasta.transdecoder.pep slirp_transcript_ids.txt > slirp_isoforms.fasta

~/git/genomeGTFtools/pfampipeline.py slirp_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001


MSTRG.10733.1
MSTRG.10733.2
ENSCGRT00000001994

grep -E 'MSTRG.10733.1|MSTRG.10733.2|ENSCGRT00000001994' stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > dnajc15_transcript_ids.txt
seqtk subseq stringtie_transcripts.fasta.transdecoder.pep dnajc15_transcript_ids.txt > dnajc15_isoforms.fasta

~/git/genomeGTFtools/pfampipeline.py dnajc15_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001

grep -E 'MSTRG.16803.1|MSTRG.16803.2|MSTRG.16803.3|ENSCGRT00000023416|ENSCGRT00000023417|ENSCGRT00000023418' stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > hikeshi_transcript_ids.txt
seqtk subseq stringtie_transcripts.fasta.transdecoder.pep hikeshi_transcript_ids.txt > hikeshi_isoforms.fasta

~/git/genomeGTFtools/pfampipeline.py hikeshi_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001


MSTRG.3523.2|MSTRG.3523.3|MSTRG.3523.4|MSTRG.3523.5|ENSCGRT00000024379
grep -E 'MSTRG.3523.2|MSTRG.3523.3|MSTRG.3523.4|MSTRG.3523.5|ENSCGRT00000024379' stringtie_transcripts.fasta.transdecoder.pep | sed 's/>//g' > immp1l_transcript_ids.txt
seqtk subseq stringtie_transcripts.fasta.transdecoder.pep immp1l_transcript_ids.txt > immp1l_isoforms.fasta

~/git/genomeGTFtools/pfampipeline.py immp1l_isoforms.fasta -P /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.hmm -c /mnt/HDD2/colin/rnaseq_manuscript_analysis/alt_splicing_analysis/pfam_db/Pfam-A.clans.tsv -e 0.0001