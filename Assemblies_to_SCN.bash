#!/bin/bash

#The output of this program is $SUBANALYSIS_$LOCUS.single.copy.framed_contigsrr.fa and is a fasta file of all the single copy genes' contigs for a particular locus
#And set of samples. This file can be aligned for gene tree building or concatenation or combined with same files from other sample sets.

#Adjust the following variables for each locus and sampleset
LOCUS=001
SUBANALYSIS=Moraceae
DIR=Test_Data/ #Name of directory where input files are
#Requirements to run:
#All de novo assembled contigs to be included in the analysis should be in a folder named Assemblies/ in DIR. There can be assemblies from other subanalyses and loci 
#In Assemblies/ as long as LOCUS is set and a file of the samples of interest is passed to -T
module load exonerate gcc/8.2.0 phyx mafft fasttree atram
#A file of sample names to be passed to -T in DIR.
#A reference sequence fasta in DIR. The name of the sequence (header) in the fasta file must match a string in the name of the assemblies

atram_framer.py -T ${DIR}${SUBANALYSIS}_Taxa.txt -a ${DIR}Assemblies/ -r ${DIR}Locus_${LOCUS}.fasta -o ${DIR}framed_contigs.${SUBANALYSIS} -l ${DIR}framer.Locus_${LOCUS}.log -m 350
#Impose a coverage cutoff to decide what sequences to drop from the framed contigs
#{1..n} where n = coverage cutoff, keeping in mind that spades header coverage is ~0.4X coverage (so n=2 is 4x coverage, roughly). change if desired.
for i in {1..19}; do grep c${i}s ${DIR}framed_contigs.${SUBANALYSIS}.Locus_${LOCUS}.fasta >> ${DIR}${SUBANALYSIS}.Drop_seqs.txt; done

#Use phyx to remove the low cov contigs
sed -i 's/>//g' ${DIR}${SUBANALYSIS}.Drop_seqs.txt
pxrms -s ${DIR}framed_contigs.${SUBANALYSIS}.Locus_${LOCUS}.fasta -f ${DIR}${SUBANALYSIS}.Drop_seqs.txt -o ${DIR}framed_contigs.${SUBANALYSIS}.Locus_${LOCUS}.fasta.rms

#Align framed contig with mafft and build tree with fasttree. Software can be changed.
mafft ${DIR}framed_contigs.${SUBANALYSIS}.Locus_${LOCUS}.fasta.rms > ${DIR}framed_contigs.Locus_${LOCUS}.fasta.rms.aln
sed -i '/^>/! s/n/-/g' ${DIR}framed_contigs.Locus_${LOCUS}.fasta.rms.aln
pxclsq -p 0.1 -s ${DIR}framed_contigs.Locus_${LOCUS}.fasta.rms.aln -o ${DIR}framed_contigs.Locus_${LOCUS}.fasta.aln-cln
FastTree -nt ${DIR}framed_contigs.Locus_${LOCUS}.fasta.aln-cln > ${DIR}framed_contigs.Locus_${LOCUS}.fasta.aln-cln.tre

module load python/2.7.14 
#(older version of python is needed for the scripts below)
#Prune tree based on branch lengths (Yang and Smith)
#change relative-cutoff and absolute-cutoff after viewing tree
python /blue/soltis/share/Greg_HR_Nitfix/phylogenomic_dataset_construction/scripts/trim_tips.py ${DIR} tre 5 .3
#Reduce clades that are all one taxon to a single tip (masking)
python /blue/soltis/share/Greg_HR_Nitfix/phylogenomic_dataset_construction/scripts/mask_tips_by_taxonID_transcripts.py ${DIR} ${DIR} y
#Write fasta file for the tips that remain in the pruned tree
python /blue/soltis/share/Greg_HR_Nitfix/phylogenomic_dataset_construction/scripts/write_fasta_files_from_trees.py ${DIR}framed_contigs.Locus_${LOCUS}.fasta.aln-cln ${DIR} mm ${DIR}

#Make csv of counts of how many contigs are left per sample
grep ">" ${DIR}framed_contigs.Locus_${LOCUS}.fasta.aln-cln|cut -d "@" -f 1|sort|uniq > ${DIR}Taxa.rr.txt
while read file; do echo -ne Locus_${LOCUS},${file}, && grep -c $file ${DIR}framed_contigsrr.fa ; done < ${DIR}Taxa.rr.txt|sed 's/>//g' > ${DIR}${SUBANALYSIS}_${LOCUS}.Taxa_contigs_per_locus_counts.csv

#Remove samples with more than one contigs from the alignment
for i in `awk -F "," '{if ($3=="1") print $2;}' ${DIR}${SUBANALYSIS}_${LOCUS}.Taxa_contigs_per_locus_counts.csv`
do grep -A1 $i ${DIR}framed_contigsrr.fa|cut -d "@" -f 1 >> ${DIR}${SUBANALYSIS}_${LOCUS}.single.copy.framed_contigsrr.fa; done
