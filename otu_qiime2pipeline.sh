#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --account=def-baum
#SBATCH --mem 50G
#SBATCH --mail-user=hepstein@uvic.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name=KI_qiime2

module load nixpkgs/16.09 qiime2/2019.4

#qiime vsearch cluster-features-closed-reference \
#  --i-table /home/hepstein/scratch/comparison/table_merged.qza \
#  --i-sequences /home/hepstein/scratch/comparison/rep-seqs_merged.qza \
#  --i-reference-sequences /home/hepstein/projects/def-baum/hepstein/KI_16S/Classifier/SILVA/ref-seqs.qza \
#  --p-perc-identity 0.97 \
#  --o-clustered-table /home/hepstein/scratch/comparison/table_97.qza \
#  --o-clustered-sequences /home/hepstein/scratch/comparison/clust_seqs_97.qza \
#  --o-unmatched-sequences /home/hepstein/scratch/comparison/unmatched.qza

#qiime feature-table summarize \
#  --i-table /home/hepstein/scratch/comparison/table_97.qza \
#  --o-visualization /home/hepstein/scratch/comparison/table_97.qzv \
#  --m-sample-metadata-file /home/hepstein/scratch/comparison/metadata_comp.txt

##Create a phylogenetic tree

#qiime phylogeny align-to-tree-mafft-fasttree \
#  --i-sequences /home/hepstein/scratch/comparison/clust_seqs_97.qza \
#  --o-alignment /home/hepstein/scratch/comparison/aligned-rep-seqs_97.qza \
#  --o-masked-alignment /home/hepstein/scratch/comparison/masked-aligned-rep-seqs_97.qza \
#  --o-tree /home/hepstein/scratch/comparison/unrooted-tree_97.qza \
#  --o-rooted-tree /home/hepstein/scratch/comparison/rooted-tree_all_97.qza

#Assign Taxonomy

qiime feature-classifier classify-sklearn \
  --i-classifier /home/hepstein/projects/def-baum/hepstein/KI_16S/Classifier/SILVA/classifier_silva.qza \
  --i-reads /home/hepstein/scratch/comparison/clust_seqs_97.qza \
  --o-classification /home/hepstein/scratch/comparison/taxonomy_97.qza

##A bug was found in the output from feature-classifier classify-sklearn where the taxon column has trailing white spaces
#e.g., "D_0__Bacteria;D-1__Firmicutes;D_2__Bacilli;D_3__Lactobacillales. " <=trailing space at end of string
#Use the following code to remove the trailing white spaces in the taxonomy artefact:

qiime tools export \
  --input-path /home/hepstein/scratch/comparison/taxonomy_97.qza \
  --output-path /home/hepstein/scratch/comparison/taxonomy-with-spaces_97
  
qiime metadata tabulate \
  --m-input-file /home/hepstein/scratch/comparison/taxonomy-with-spaces_97/taxonomy.tsv  \
  --o-visualization /home/hepstein/scratch/comparison/taxonomy-with-spaces_97/taxonomy-as-metadata.qzv

qiime tools export \
 --input-path /home/hepstein/scratch/comparison/taxonomy-with-spaces_97/taxonomy-as-metadata.qzv \
 --output-path /home/hepstein/scratch/comparison/taxonomy-with-spaces_97/taxonomy-as-metadata

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-path /home/hepstein/scratch/comparison/taxonomy-with-spaces_97/taxonomy-as-metadata/metadata.tsv \
  --output-path /home/hepstein/scratch/comparison/tax_97.qza

#Filter out mitochondria, chloroplasts, eukaryotes
#Filter table 

qiime taxa filter-table \
  --i-table /home/hepstein/scratch/comparison/table_97.qza \
  --i-taxonomy /home/hepstein/scratch/comparison/tax_97.qza \
  --p-exclude chloroplast,eukaryote \
 --o-filtered-table /home/hepstein/scratch/comparison/table_97_filtered.qza
   
qiime feature-table summarize \
  --i-table /home/hepstein/scratch/comparison/table_97_filtered.qza \
  --o-visualization /home/hepstein/scratch/comparison/table_97_filtered.qzv

#Check for alpha rarefaction at 5000 reads
qiime diversity alpha-rarefaction \
  --i-table /home/hepstein/scratch/comparison/table_97_filtered.qza \
  --i-phylogeny /home/hepstein/scratch/comparison/rooted-tree_all_97.qza \
  --p-max-depth 5000 \
  --m-metadata-file /home/hepstein/scratch/comparison/metadata_comp.txt \
  --o-visualization /home/hepstein/scratch/comparison/97_alpha-rarefaction.qzv

##Make taxonomic barplots to view in view.qiime

qiime taxa barplot \
  --i-table /home/hepstein/scratch/comparison/table_97_filtered.qza \
  --i-taxonomy /home/hepstein/scratch/comparison/tax_97.qza \
  --m-metadata-file /home/hepstein/scratch/comparison/metadata_comp.txt \
  --o-visualization /home/hepstein/scratch/comparison/97_barplot.qzv
  
  





