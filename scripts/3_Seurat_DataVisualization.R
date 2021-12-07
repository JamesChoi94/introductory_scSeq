
###### Data Visualization using Seurat functions ######

# Load libraries
library('Seurat')
library('ggplot2')
library('dplyr')

macroglia <- readRDS(file = './data/macroglia.rds')



# Seurat Data Organization ------------------------------------------------


# Gene expression data are preprocessed/normalized before they are clustered. 
# The raw data, normalized data, scaled data, and meta-data are all stored in
# separate "slots" of the Seurat data structure.

# Use the "@" symbol to access the different slots.


# Gene expression values are stored in the "assays" slot. Under the "assays"
# slot, there are 3 assays. Please use the "RNA" assay. The "SCT" and 
# "integrated" assays are used internally for batch correction and should
# not be used for visualization/differential gene expression tests. To use the
# "RNA" assay requires setting the DefaultAssay(). This is a one-time operation.
# If further analysis requiring the "SCT" or "integrated" is needed, you will
# need to change the DefaultAssay before doing so. 

DefaultAssay(macroglia) <- 'RNA'

macroglia@assays$RNA
macroglia@assays$RNA@counts # raw mRNA count values
macroglia@assays$RNA@data # Normalized for total mRNA content per cell
macroglia@assays$RNA@scale.data # Relative expression compared to all other cells (scaled by standard deviation)


# Relative difference vs absolute difference
# 1, 2, 3
#    vs
# 10, 20, 30


dim(macroglia@assays$RNA@counts)
macroglia@assays$RNA@counts['Gfap',]
summary(macroglia@assays$RNA@counts['Gfap',])


# Important meta data (characteristics of individual cells e.g. Time) can 
# be found under the "meta.data" slot. This data is stored in a data.frame, 
# where different characteristics (e.g. Time, cell-type, etc.) run along
# the columns and the values for each individual cell run along the rows.

macroglia@meta.data
dim(macroglia@meta.data)
colnames(macroglia@meta.data) # what categorical variables can we choose from?
rownames(macroglia@meta.data)


# The results of the subcluster analysis (as reported in our pre-print) are
# saved under the "subcluster" column of the macroglia@meta.data data.

macroglia@meta.data$macroglia_subcluster
table(macroglia@meta.data$macroglia_subcluster)
table(macroglia@meta.data$macroglia_subcluster, macroglia@meta.data$time)


# There are other numerical/categorical data that you can explore.

table(macroglia@meta.data$Phase)
summary(macroglia@meta.data$percent.mt)


# At any one time, Seurat assigns each cell an "identity". These are the 
# labels that are currently set to each cell. For example:

Idents(macroglia)
table(Idents(macroglia))


# You can set different cell labels using the data stored in the meta.data.

macroglia@meta.data$time
Idents(macroglia) <- 'time' # Write meta.data column name here
Idents(macroglia)
table(Idents(macroglia))

# Taking a subset of cells:
my_subset <- macroglia[,macroglia@meta.data$time == '1dpi']
dim(my_subset)



# UMAP Visualizations -----------------------------------------------------

# Seurat has many visualization functions. See: https://satijalab.org/seurat/v3.1/visualization_vignette.html

# "DimPlot" stands for Dimensional (reduction) plot. 
# To generate a UMAP plot:
DimPlot(macroglia)


# To overlay the cell identities (labels):
Idents(macroglia) <- 'macroglia_subcluster'
DimPlot(macroglia, label = TRUE)
DimPlot(macroglia, label = TRUE, label.size = 5)


# To remove the legend on the side:
DimPlot(macroglia) + NoLegend()
DimPlot(macroglia, label = TRUE, label.size = 5) + NoLegend()

# To increase dot size:
DimPlot(macroglia, pt.size = 2)

# You can change how the cells are colored by setting the identity.
Idents(macroglia) <- 'Time'
DimPlot(macroglia)

Idents(macroglia) <- 'orig.ident' # sample of origin, we had duplicates per time
DimPlot(macroglia)


# You can also generate separate UMAPs side-by-side.
# Generate a separate UMAP for each time-point:
Idents(macroglia) <- 'macroglia_subcluster'
DimPlot(macroglia, split.by = 'time')
DimPlot(macroglia, split.by = 'orig.ident', ncol = 4)
DimPlot(macroglia, group.by = 'macroglia_subcluster', split.by = 'time')



# You can store the output plots as a variable:
Idents(macroglia) <- 'subcluster'
plot_by_time <- DimPlot(macroglia, label = TRUE, label.size = 5) + NoLegend()
plot_by_time


# Two ways to save a plot:
# 1) Click "Export" in the plots pane.
# 2) Run the following line:
ggsave(plot = plot_by_time, 
       filename = 'macroglia_by_subcluster.tiff', 
       device = 'tiff',
       height = 5.5, width = 6.5)



# The above plots had cells colored by a categorial variable. How about for
# a numerical variable? Use the FeaturePlot() function.

# Plot gene expression values:
FeaturePlot(macroglia, features = 'Aqp4')
FeaturePlot(macroglia, features = 'Aqp4', split.by = 'time')
FeaturePlot(macroglia, features = 'Aqp4', slot = 'counts')


# Plot numerical variables from the meta.data
FeaturePlot(macroglia, 'nFeature_RNA')
FeaturePlot(macroglia, 'nFeature_RNA', order = TRUE)
FeaturePlot(macroglia, 'Slc1a2', order = TRUE)


# Plot multiple features:
FeaturePlot(macroglia, features = c('Aqp4', 'Foxj1'))
my_genes <- c('Aqp4','Foxj1','Slc1a2','Gfap','Olig2')
FeaturePlot(macroglia, features = my_genes)


# Change colors if desired:
FeaturePlot(macroglia, features = 'Foxj1', cols = c('purple','yellow'))


# Be sure to write the correct gene name:
FeaturePlot(macroglia, 'Cd11b')
FeaturePlot(macroglia, 'Itgam')



# More quantitative visualization of gene expression ----------------------


# Depending on your goal, you may want to use a different plot. For example,
# if we want to show relative expression of a gene, we may consider a violin
# plot:

VlnPlot(macroglia, 'Aqp4')
VlnPlot(macroglia, 'Aqp4', split.by = 'time')

# color violins manually
my_violin_colors <- c('Ependymal-A' = 'red',
                      'Ependymal-B' = 'blue',
                      'Astroependymal' = 'green',
                      'Astrocyte' = '#fff000',
                      'OPC-A' = 'magenta',
                      'OPC-B' = 'navy',
                      'Pre-Oligo' = 'purple',
                      'Oligodendrocyte' = 'grey')
VlnPlot(macroglia, features = 'Aqp4', pt.size = 0) +
  scale_fill_manual(values = my_violin_colors)

# Plot specific subtypes and split by condition
Idents(macroglia) <- 'macroglia_subcluster'
my_favorite_cells <- c('Astrocyte','OPC-A','OPC-B')
VlnPlot(macroglia, features = 'Apoe', split.by = 'time',
        idents = my_favorite_cells) +
  scale_fill_manual(values = c('Uninjured' = 'red', 
                               '1dpi' = 'yellow', 
                               '3dpi' = 'green', 
                               '7dpi' = 'blue'))


Idents(macroglia) <- 'time'
VlnPlot(macroglia, 'Apoe')
VlnPlot(macroglia, 'Apoe', split.by = 'subcluster')


Idents(macroglia) <- 'subcluster'
VlnPlot(macroglia, features = 'Apoe', split.by = 'orig.ident', idents = c('Astrocyte', 'OPC'))
VlnPlot(macroglia, 'Apoe', pt.size = 0)



# Scaled Expression -------------------------------------------------------

# Scaled expression is a z-scored relative expression value. It considers
# the variance of the gene's normalized expression value and scales each
# cell's expression such that the standard deviation is 1.

# For visualizations where a several gene's are being mapped to the same
# scale, expressions need to be scaled.


# You can manually enter any genes you would like to scale:
my_genes <- c('Gfap','Aqp4','Rgs5','Mki67')
macroglia <- ScaleData(macroglia, features = my_genes)

# Expression values will now take negative numbers.
macroglia@assays$RNA@scale.data['Mki67',]
summary(macroglia@assays$RNA@scale.data['Mki67',])




# Visualizing relative expression -----------------------------------------

# First, we can check via dot plots.
DotPlot(macroglia, features = my_genes)


# We can further split the dots by other meta.data (e.g. Time)
DotPlot(macroglia, 
        features = my_genes,
        split.by = 'time',
        cols = c('red','red','red','red'))


# Rotate axis for better spacing
DotPlot(macroglia, 
        features = c('Slc1a2','Aqp4','Apoe'),
        split.by = 'time',
        cols = rep('red', 4)) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(macroglia, features = my_genes) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# We can further inspect relative expression at the single-cell levels with
# a heatmap:
DoHeatmap(macroglia, features =  my_genes)

DoHeatmap(macroglia, features =  my_genes, size = 3)



# Save Data ---------------------------------------------------------------

# If there are any changes we have made that we would like to save/store in the
# Seurat object, we will need to save the R-readable dataset as an .rds file.

# WARNING: this will overwrite the previous "macroglia.rds" file.
saveRDS(macroglia, 'macroglia.rds')
