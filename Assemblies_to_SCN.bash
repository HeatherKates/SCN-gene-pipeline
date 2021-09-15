#!/bin/bash

#This file should be in each subanalysis/locus dir

#Adjust the following variables for each locus and sampleset
LOCUS=
SUBANALYSIS=

#Requirements to run:
#All de novo assembled contigs to be included in the analysis should be in a folder named Assemblies/ in this dir
#module load python3 exonerate gcc/8.2.0 phyx mafft fasttree atram
#aTRAM/atram_framer.py
#A file of sample names to be passed to -T in ../
#A reference sequence fasta in this dir

atram_framer.py -T ../${SUBANALYSIS}_Taxa.txt -a Assemblies/ -r Locus_${LOCUS}.fasta -o framed_contigs.${SUBANALYSIS} -l framer.Locus_${LOCUS}.log -m 350
#Impose a coverage cutoff to decide what sequences to drop from the framed contigs
#{1..n} where n = coverage cutoff, keeping in mind that spades header coverage is ~0.4X coverage (so n=2 is 4x coverage, roughly). change if desired.
for i in {1..19}; do grep c${i}s framed_contigs.${SUBANALYSIS}.Locus_${LOCUS}.fasta >> ${SUBANALYSIS}.Drop_seqs.txt; done

#Use phyx to remove the low cov contigs
sed -i 's/>//g' ${SUBANALYSIS}.Drop_seqs.txt
pxrms -s framed_contigs.${SUBANALYSIS}.Locus_${LOCUS}.fasta -f ${SUBANALYSIS}.Drop_seqs.txt -o framed_contigs.${SUBANALYSIS}.Locus_${LOCUS}.fasta.rms

#Align framed contig with mafft and build tree with fasttree. Software can be changed.
mafft framed_contigs.Locus_${LOCUS}.fasta.rms > framed_contigs.Locus_${LOCUS}.fasta.rms.aln
sed -i '/^>/! s/n/-/g' framed_contigs.Locus_${LOCUS}.fasta.rms.aln
pxclsq -p 0.1 -s framed_contigs.Locus_${LOCUS}.fasta.rms.aln -o framed_contigs.Locus_${LOCUS}.fasta.aln-cln
FastTree -nt framed_contigs.Locus_${LOCUS}.fasta.aln-cln > framed_contigs.Locus_${LOCUS}.fasta.aln-cln.tre

#module load python/2.7.14 (older version of python is needed for the scripts below)
#Prune tree based on branch lengths (Yang and Smith)
#change relative-cutoff and absolute-cutoff after viewing tree
python /ufrc/soltis/share/Greg_HR_Nitfix/phylogenomic_dataset_construction/scripts/trim_tips.py . tre 5 .3
#Reduce clades that are all one taxon to a single tip (masking)
python /ufrc/soltis/share/Greg_HR_Nitfix/phylogenomic_dataset_construction/scripts/mask_tips_by_taxonID_transcripts.py . . y
#Write fasta file for the tips that remain in the pruned tree
python /ufrc/soltis/share/Greg_HR_Nitfix/phylogenomic_dataset_construction/scripts/write_fasta_files_from_trees.py framed_contigs.Locus_${LOCUS}.fasta.aln-cln . mm .

#Make csv of counts of how many contigs are left per sample
grep ">" framed_contigs.Locus_${LOCUS}.fasta.aln-cln|cut -d "@" -f 1|sort|uniq > Taxa.rr.txt
while read file; do echo -ne Locus_${LOCUS},${file}, && grep -c $file framed_contigsrr.fa ; done < Taxa.rr.txt|sed 's/>//g' > ${SUBANALYSIS}_${LOCUS}.Taxa_contigs_per_locus_counts.csv

#Remove samples with more than one contigs from the alignment
for i in `awk -F "," '{if ($3=="1") print $2;}' ${SUBANALYSIS}_${LOCUS}.Taxa_contigs_per_locus_counts.csv`
do grep -A1 $i framed_contigsrr.fa|cut -d "@" -f 1 >> ${SUBANALYSIS}_${LOCUS}.single.copy.framed_contigsrr.fa; done
