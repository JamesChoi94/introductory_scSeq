
# For this demo, we'll be working with a dataset published by Li et al from Zhigang He's laboratory. Dataset contains single-cell RNAseq data from FACS-sorted microglia. 

# GEO accession: GSE150871

# Link to GEO submission: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150871

# Link to paper: https://doi.org/10.1038/s41586-020-2795-6

# We will be following the Seurat pipeline as outlined in their tutorial here: https://satijalab.org/seurat/articles/integration_introduction.html


# Libraries and directories -----------------------------------------------


# Check working directory. Change if needed.
getwd()


# Load libraries
require('Seurat')
require('ggplot2')
require('dplyr')


# Data download -----------------------------------------------------------


# Download directly from my Google drive: https://drive.google.com/drive/folders/1LpFPA-a_kaqGevxqF5jFldptOF8BrkmK?usp=sharing

# The .rds files contain the count matrices for each of the 5 samples from the study. Rows are genes and columns are cells. They are saved as dgCmatrix objects. 

counts_files <- list.files(path = 'D:/MiamiProject/Microglia_Fate/data/counts/',
                           full.names = TRUE)
short_file_names <- list.files(
  path = 'D:/MiamiProject/Microglia_Fate/data/counts/',
  full.names = FALSE
)
counts <- vector(mode = 'list', length = length(counts_files))
for (i in 1:length(counts_files)) {
  counts[[i]] <- readRDS(file = counts_files[i])
}


# Quality Control ---------------------------------------------------------


# Calculate percent.mt and percent.rp for each data set
for (i in 1:length(counts)) {
  counts[[i]] <- CreateSeuratObject(
    counts = counts[[i]],
    project = short_file_names[i]
  )
  counts[[i]] <- PercentageFeatureSet(
    object = counts[[i]],
    pattern = '^mt-',
    col.name = 'percent.mt'
  )
  counts[[i]] <- PercentageFeatureSet(
    object = counts[[i]],
    pattern = '^Rp[ls]',
    col.name = 'percent.rp'
  )
}

# Extract these fields from Seurat object metadata
meta_feats <- c('orig.ident','nCount_RNA','nFeature_RNA','percent.mt','percent.rp')
metadata <- list()
for (i in 1:length(counts)) {
  metadata[[i]] <- counts[[i]]@meta.data[meta_feats]
}

# Merge metadata data.frames into a single data.frame (this helps automatically set common scale bars for plots)
metadata <- Reduce(rbind, metadata)
p1 <- metadata %>%
  reshape2::melt(id.vars = c('orig.ident')) %>%
  ggplot(mapping = aes(x = orig.ident, y = value)) + 
  geom_violin(mapping = aes(fill = orig.ident), scale = 'width') +
  facet_wrap(. ~ variable, scales = 'free_y') +
  theme_bw() +
  xlab(label = 'Sample ID') +
  ylab(label = 'Value') +
  theme(axis.text.x = element_blank())
p1

# Calculate median-absolute-deviation based thresholds for nCount_RNA, which allows for variability in mean sequencing depth by sample which can be affected by number of cells. Set flat threshold at 15% for percent.mt based on similarity in distribution between samples. Filter cells.  
get_mad_max <- function(x, dev = 3) {
  return(median(x) + dev * mad(x, center = 1))
}

high_qc <- vector(mode = 'list', length = length(counts))
for (i in 1:length(counts)) {
  # compute UMI count corresponding to 3 MADs above median
  umi_mad_max <- get_mad_max(counts[[i]]$nCount_RNA, dev = 3)
  high_qc[[i]] <- 
    # remove cells with UMI count above 3 MADs
    (counts[[i]]@meta.data$nCount_RNA <= umi_mad_max) &
    # remove cells with percent.mt > 15
    (counts[[i]]@meta.data$percent.mt <= 15)
}

# Summary table
filter_results <- cbind(
  'Before' = sapply(counts, ncol),
  'Filtered' = sapply(high_qc, sum)
)
rownames(filter_results) <- short_file_names
filter_results

# Filter out cells and extract count matrix
keep_counts <- list()
for (i in 1:length(counts)) {
  keep_counts[[i]] <- counts[[i]]@assays$RNA@counts[, high_qc[[i]]]
}

# To clear RAM if needed:
# rm(counts, high_qc)


# Merged Cluster Analysis -------------------------------------------------

# Ideally, samples generated from the same study do not require batch correction or other normalization techniques to compare samples or experimental conditions. However, sometimes designing experiments in such a way to completely avoid batch effects is not practical. 

# Initially, we test whether the data requires batch correction. We do this by merging the count matrices together and running it through the Seurat pipeline. Our metric for assessing batch effects is a qualitative observation of the resulting UMAP. More quantitative approaches exist and are routinely used.

# In order to merge count matrices, we need to make sure each cell (column) has a unique identifier i.e. barcode. Since oligonucleotide barcode sequences are re-used between 10X runs, we append the sample ID to the oligonucleotide to make them unique. The sample IDs are based on the order of the count matrix file names.

short_file_names
sample_id <- c('Ctl_1', '3dpi_1', '3dpi_2', '5dpi_1', '5dpi_2')
for (i in 1:length(keep_counts)) {
  colnames(keep_counts[[i]]) <- paste(
    colnames(keep_counts[[i]]),
    sample_id[i],
    sep = '_'
  )
}

# Run through Seurat pipeline
merged_counts <- Reduce(f = cbind, x = keep_counts)
merged_counts <- CreateSeuratObject(counts = merged_counts)
merged_counts <- NormalizeData(merged_counts)
merged_counts <- FindVariableFeatures(merged_counts, nfeatures = 2500)
merged_counts <- ScaleData(merged_counts)
merged_counts <- RunPCA(merged_counts, npcs = 30)
merged_counts <- FindNeighbors(merged_counts, dims = 1:15)
merged_counts <- RunUMAP(merged_counts, dims = 1:15)
merged_counts <- FindClusters(merged_counts)

# Need to create meta.data column containing sample ID. We can extract them from the modified oligonucleotide cell barcodes from earlier.

cell_sample_id <- colnames(merged_counts)
cell_sample_id <- strsplit(cell_sample_id, split = '_')
cell_sample_id <- sapply(
  X = cell_sample_id,
  FUN = function(x) paste(x[2], x[3], sep = '_')
)
merged_counts$sample_id <- cell_sample_id
table(merged_counts$sample_id)

# Inspect UMAP of cells grouped by sample
DimPlot(merged_counts, group.by = 'sample_id', shuffle = TRUE, pt.size = 1)


# The UMAP shows that the two replicates for 3dpi overlap well and that the two replicates for 5dpi also overlap, indicating that batch effects between these technical replicates are minimal. However, we see there is a systematic shift between the Ctl sample, 3dpi samples, and 5dpi samples in the larger group of cells (which are microglia - you may run FeaturePlot(merged_counts, 'P2ry12') to check). Whether these shifts are biologically driven or technically driven is impossible to determine due to the experimental design (each sample was sequenced separately). 

# Therefore, we assume some of these shifts between time-points are due to batch effect and proceed with the Seurat Integration pipeline. Since there is no systematic shift between replicates within the time-points, these datasets can be merged without any correction. 

# To clear RAM if needed:
# rm(merged_counts)

corrected_counts <- list(
  'Ctl' = keep_counts[[1]],
  'D3' = cbind(keep_counts[[2]], keep_counts[[3]]),
  'D5' = cbind(keep_counts[[4]], keep_counts[[5]])
)
for (i in 1:length(corrected_counts)) {
  corrected_counts[[i]] <- CreateSeuratObject(
    counts = corrected_counts[[i]],
    project = names(corrected_counts)[i]
  )
  corrected_counts[[i]] <- NormalizeData(corrected_counts[[i]])
  corrected_counts[[i]] <- FindVariableFeatures(corrected_counts[[i]])
}

features <- SelectIntegrationFeatures(corrected_counts, nfeatures = 2000)
anchors <- FindIntegrationAnchors(corrected_counts, anchor.features = features)
corrected_counts <- IntegrateData(anchors)
DefaultAssay(corrected_counts) <- 'integrated'
corrected_counts <- ScaleData(corrected_counts)
corrected_counts <- RunPCA(corrected_counts)
corrected_counts <- FindNeighbors(corrected_counts, dims = 1:10)
corrected_counts <- RunUMAP(corrected_counts, dims = 1:10)
corrected_counts <- FindClusters(corrected_counts)

# re-assign sample ID as before
cell_sample_id <- colnames(corrected_counts)
cell_sample_id <- strsplit(cell_sample_id, split = '_')
cell_sample_id <- sapply(
  X = cell_sample_id,
  FUN = function(x) paste(x[2], x[3], sep = '_')
)
corrected_counts$sample_id <- cell_sample_id

# Inspect UMAP of cells grouped by sample
DimPlot(corrected_counts, group.by = 'sample_id', shuffle = TRUE)


# The resulting UMAP shows that many of the cells from the 3 time-points overlap. Again, whether this overlap is biologically driven or technically driven is currently impossible to precisly determine. However, we can inspect the UMAP split by time-point to assess whether any cells or clusters of cells are unique to any of the time-points.
DimPlot(corrected_counts, group.by = 'sample_id', split.by = 'orig.ident',
        shuffle = TRUE)

# Here, the UMAP shows slight differences in cell distribution among the larger group of cells (microglia). 
DimPlot(corrected_counts, label = TRUE, label.size = 5, shuffle = TRUE)

# Cells in cluster 9 and 14 are enriched in Ctl and present in lower numbers in the 3dpi and even lower numbers in 5dpi. The markers that define these clusters are related to cell division.
DefaultAssay(corrected_counts) <- 'RNA'
FeaturePlot(corrected_counts, 'Mki67', split.by = 'orig.ident', order = TRUE)
