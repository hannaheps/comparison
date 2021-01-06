#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --account=def-baum
#SBATCH --mem 50G
#SBATCH --mail-user=hepstein@uvic.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name=KI_qiime2

module load nixpkgs/16.09 qiime2/2019.4


#qiime cutadapt trim-single \
#  --i-demultiplexed-sequences /home/hepstein/scratch/comparison/jmi_demux.qza \
#  --p-front GTGYCAGCMGCCGCGGTAA \
#  --o-trimmed-sequences /home/hepstein/scratch/comparison/jmi_demux_trim.qza

qiime dada2 denoise-single \
  --i-demultiplexed-seqs /home/hepstein/scratch/comparison/EMP_comp_demux.qza \
  --p-trunc-len 120 \
  --o-table /home/hepstein/scratch/comparison/emp_table.qza \
  --o-representative-sequences /home/hepstein/scratch/comparison/emp_rep-seqs.qza \
  --o-denoising-stats /home/hepstein/scratch/comparison/emp_denoising-stats.qza

#qiime dada2 denoise-single \
#  --i-demultiplexed-seqs /home/hepstein/scratch/comparison/jmi_demux_trim.qza \
#  --p-trunc-len 120 \
#  --o-table /home/hepstein/scratch/comparison/table_jmi.qza \
#  --o-representative-sequences /home/hepstein/scratch/comparison/rep-seqs_jmi.qza \
#  --o-denoising-stats /home/hepstein/scratch/comparison/denoising-stats_jmi.qza
  
 qiime feature-table merge \
  --i-tables /home/hepstein/scratch/comparison/emp_table.qza \
  --i-tables /home/hepstein/scratch/comparison/table_jmi.qza \
  --o-merged-table /home/hepstein/scratch/comparison/table_merged.qza

qiime feature-table merge-seqs \
  --i-data /home/hepstein/scratch/comparison/emp_rep-seqs.qza \
  --i-data /home/hepstein/scratch/comparison/rep-seqs_jmi.qza \
  --o-merged-data /home/hepstein/scratch/comparison/rep-seqs_merged.qza
  
  
##Once merged, from here on we can use "table_all.qza" and "rep-seqs_all.qza"
##Metadata has to be in tab separated formatting

qiime feature-table summarize \
  --i-table /home/hepstein/scratch/comparison/table_merged.qza  \
  --o-visualization /home/hepstein/scratch/comparison/table_merged.qza  \
  --m-sample-metadata-file /home/hepstein/scratch/comparison/metadata_comp.txt

qiime feature-table tabulate-seqs \
  --i-data /home/hepstein/scratch/comparison/rep-seqs_merged.qza  \
  --o-visualization /home/hepstein/scratch/comparison/rep-seqs_merged.qzv


#Use the classifier to assign taxonomy - this can take a long time and needs lots of memory (50G?) 

qiime feature-classifier classify-sklearn \
  --i-classifier /home/hepstein/projects/def-baum/hepstein/KI_16S/Classifier/SILVA/classifier_silva.qza \
  --i-reads /home/hepstein/scratch/comparison/rep-seqs_merged.qza \
  --o-classification /home/hepstein/scratch/comparison/tax.qza

##A bug was found in the output from feature-classifier classify-sklearn where the taxon column has trailing white spaces
#e.g., "D_0__Bacteria;D-1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales. " <=trailing space at end of string
#Use the following code to remove the trailing white spaces in the taxonomy artefact:

qiime tools export \
  --input-path /home/hepstein/scratch/comparison/tax.qza \
  --output-path /home/hepstein/scratch/comparison/taxonomy-with-spaces
  
qiime metadata tabulate \
  --m-input-file /home/hepstein/scratch/comparison/taxonomy-with-spaces/taxonomy.tsv  \
  --o-visualization /home/hepstein/scratch/comparison/taxonomy-with-spaces/taxonomy-as-metadata.qzv

qiime tools export \
 --input-path /home/hepstein/scratch/comparison/taxonomy-with-spaces/taxonomy-as-metadata.qzv \
 --output-path /home/hepstein/scratch/comparison/taxonomy-as-metadata

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path /home/hepstein/scratch/comparison/taxonomy-as-metadata/metadata.tsv \
  --output-path /home/hepstein/scratch/comparison/tax-without-spaces.qza

#The following code spits out a table of features, or ASVs, present in sample set
#You can view by dragging the .qzv file into view.qiime2.org
#Or use the "qiime tools view" command

qiime metadata tabulate \
  --m-input-file /home/hepstein/scratch/comparison/tax-without-spaces.qza \
  --o-visualization /home/hepstein/scratch/comparison/tax-without-spaces.qzv
  
### Filtering your table: removing mitochondria and chloroplast reads, PCR/sequencing errors, etc###

##Filter table

qiime taxa filter-table \
  --i-table /home/hepstein/scratch/comparison/table_merged.qza \
  --i-taxonomy /home/hepstein/scratch/comparison/tax-without-spaces.qza \
  --p-exclude chloroplast,eukaryote \
  --o-filtered-table /home/hepstein/scratch/comparison/table_merged_filtered.qza
   
qiime feature-table summarize \
  --i-table /home/hepstein/scratch/comparison/table_merged_filtered.qza \
  --o-visualization /home/hepstein/scratch/comparison/table_merged_filtered.qzv
  
##Create a phylogenetic tree (With large datasets this takes ALOT of RAM - give it 50G)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /home/hepstein/scratch/comparison/rep-seqs_merged.qza \
  --o-alignment /home/hepstein/scratch/comparison/aligned-rep-seqs.qza \
  --o-masked-alignment /home/hepstein/scratch/comparison/masked-aligned-rep-seqs.qza \
  --o-tree /home/hepstein/scratch/comparison/unrooted-tree.qza \
  --o-rooted-tree /home/hepstein/scratch/comparison/rooted-tree.qza
 
  