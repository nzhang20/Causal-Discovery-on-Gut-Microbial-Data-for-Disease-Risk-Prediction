#!/bin/bash

# Activate the QIIME2 environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-2024.10

# Set input and output paths
MANIFEST_PATH="manifest.tsv"
METADATA_PATH="our-metadata.tsv"

# Import data commands for single-end with manifest
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path "$MANIFEST_PATH" \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

# Use this visual to find your denoising parameters
qiime demux summarize \
  --i-data single-end-demux.qza \
  --o-visualization single-end-demux.qzv

# Denoising 1: Quality score filtering
qiime quality-filter q-score \
 --i-demux single-end-demux.qza \
 --o-filtered-sequences demux-filtered.qza \
 --o-filter-stats demux-filter-stats.qza

qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv

# Denoising 2, with Deblur:
# adjust trim length + jobs/cores to be used
qiime deblur denoise-16S \
 --i-demultiplexed-seqs demux-filtered.qza \
 --p-trim-length 150 \
 --p-jobs-to-start 12 \
 --o-representative-sequences rep-seqs-150.qza \
 --o-table table-150.qza \
 --p-sample-stats \
 --o-stats deblur-stats-150.qza

# Generate artifacts for visualization
qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats-150.qza \
  --o-visualization deblur-stats-150.qzv

qiime feature-table summarize \
  --i-table table-150.qza \
  --o-visualization table-150.qzv \
  --m-sample-metadata-file "$METADATA_PATH"

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-150.qza \
  --o-visualization rep-seqs-150.qzv

#Start tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-150.qza \
  --p-n-threads 12 \
  --o-alignment aligned-rep-seqs-150.qza \
  --o-masked-alignment masked-aligned-rep-seqs-150.qza \
  --o-tree unrooted-tree-150.qza \
  --o-rooted-tree rooted-tree-150.qza

# Select sample-depth
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree-150.qza \
  --i-table table-150.qza \
  --p-n-jobs-or-threads 12 \
  --p-sampling-depth 2000 \
  --m-metadata-file "$METADATA_PATH" \
  --output-dir core-metrics-results-2000-150bp

# Core metrics visualization

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-2000-150bp/faith_pd_vector.qza \
  --m-metadata-file "$METADATA_PATH" \
  --o-visualization core-metrics-results-2000-150bp/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-2000-150bp/evenness_vector.qza \
  --m-metadata-file "$METADATA_PATH" \
  --o-visualization core-metrics-results-2000-150bp/evenness-group-significance.qzv


qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-2000-150bp/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file "$METADATA_PATH" \
  --m-metadata-column t2d \
  --o-visualization core-metrics-results-2000-150bp/unweighted-unifrac-t2d-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results-2000-150bp/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file "$METADATA_PATH" \
  --m-metadata-column country \
  --o-visualization core-metrics-results-2000-150bp/unweighted-unifrac-country-group-significance.qzv \
  --p-pairwise

#Emperor plots, if continuous parameter included:
qiime emperor plot \
  --i-pcoa core-metrics-results-2000-150bp/bray_curtis_pcoa_results.qza \
  --m-metadata-file "$METADATA_PATH" \
  --p-custom-axes days-since-experiment-start \
  --o-visualization core-metrics-results-2000-150bp/bray-curtis-emperor-days-since-experiment-start.qzv


#Alpha rarefaction plotting
qiime diversity alpha-rarefaction \
  --i-table table-150.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 4000 \
  --m-metadata-file "$METADATA_PATH" \
  --o-visualization alpha-rarefaction-4000.qzv

# Taxonomy
qiime feature-classifier classify-sklearn \
  --i-classifier  2024.09.backbone.full-length.nb.sklearn-1.4.2.qza \
  --i-reads rep-seqs-150.qza \
  --p-n-jobs 0 \
  --o-classification downstream/taxonomy-gg-2024-09-full.qza

qiime metadata tabulate \
  --m-input-file downstream/taxonomy-gg-2024-09-full.qza \
  --o-visualization downstream/taxonomy-gg-2024-09-full.qzv

qiime taxa barplot \
  --i-table table-150.qza \
  --i-taxonomy downstream/taxonomy-gg-2024-09-full.qza \
  --m-metadata-file "$METADATA_PATH" \
  --o-visualization downstream/taxa-bar-plots-gg-2024-09.qzv