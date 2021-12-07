

##### Differential Gene Expression (DE) tests using Seurat #####

# Libraries
library('Seurat')
library('ggplot2')
library('dplyr')
install.packages('ggrepel') # if you do not have 'ggrepel' installed.
library('ggrepel')

# Load data
macroglia <- readRDS(file = 'data/macroglia.rds')

# DE tests _must_ be done with the RNA assay. 
DefaultAssay(macroglia) <- 'RNA'

# DE tests compare between identities
Idents(macroglia) <- 'subcluster'



# Calculating Differential Expression -------------------------------------


# FindMarkers() is the basic DE test function. The output is a data.frame that 
# contains results for all of the genes that were tested. The column data are:
# 
#   1. p_val: p-values for Wilcoxon Rank Sum test (default test)
#   2. avg_logFC: average log2(fold-change). Fold-change calculated as the ratio
#     between average expression in first group over average expression in second
#     group. Positive values indicate higher expression in first group (ident.1).
#   3. pct.1: Percent of cells in ident.1 that express a gene.
#   4. pct.2: Percent of cells in ident.2 that express a gene.
#   5. p_val_adj: Adjusted p-value for multiple hypothesis correction. Uses 
#     bonferroni correction (default).
#     
# For more detailed info, run "?FindMarkers()"


# To calculate differentially expressed genes for a cluster vs all other cells:
my_de <- FindMarkers(macroglia,
                     ident.1 = 'Astroependymal')


# To compare two cell-types:
my_de <- FindMarkers(macroglia,
                     ident.1 = 'Astroependymal', # required
                     ident.2 = 'Astrocyte')


# To compare each cluster with all other clusters for all clusters (to identify
# marker genes, usually we desire positive markers, so we set "only.pos" to TRUE):
my_de <- FindAllMarkers(macroglia, only.pos = TRUE)


# To compare one cell-type with multiple cell-types (it will combine all cells 
# as one group):
my_de <- FindMarkers(macroglia,
                     ident.1 = 'Astroependymal',
                     ident.2 = c('Ependymal-A','Ependymal-B','Ependymal-C'))


# By default, only genes that have a minimum log(fold-change) of at least 0.25 
# between the two groups are tested. This is partially to reduce computation 
# time and partially to reduce false positives. To adjust this threshold, set the
# "logfc.threshold" argument:
my_de <- FindMarkers(macroglia,
                     ident.1 = 'Astroependymal',
                     ident.2 = 'Astrocyte',
                     logfc.threshold = 0)


# By default, only genes that have a minimum percent expression of 10% in either
# of the two groups are tested. This is partially to reduce computation time and
# partially to reduce false positives. To adjust this threshold, set the
# "min.pct" argument:
my_de <- FindMarkers(macroglia,
                     ident.1 = 'Astroependymal',
                     ident.2 = 'Astrocyte',
                     min.pct = 0.2)

# Test all ~33000 genes in dataset (whole reference genome, mm10)
my_de <- FindMarkers(macroglia,
                     ident.1 = 'Astroependymal',
                     ident.2 = 'Astrocyte',
                     logfc.threshold = 0,
                     min.pct = 0)


# You can also split a certain cell-type by Time and perform DE between time 
# points.
#   1) Specify the cell-type(s) with "subset.ident"
#   2) Set "group.by" to "Time"
#   3) Select the two groups (time-points)
my_de <- FindMarkers(macroglia,
                     subset.ident = 'Astrocyte',
                     group.by = 'Time',
                     ident.1 = '3d',
                     ident.2 = '1d')
my_de <- FindMarkers(macroglia,
                     subset.ident = c('Ependymal-A','Ependymal-B','Ependymal-C'),
                     group.by = 'Time',
                     ident.1 = '7d',
                     ident.2 = '3d')


# Conversely, take all cells from a time point, then group cells by subcluster
# identities, finally compare two groups.
Idents(macroglia) <- 'Time'
my_de <- FindMarkers(macroglia,
                     subset.ident = '1d',
                     group.by = 'subcluster',
                     ident.1 = 'Astroependymal',
                     ident.2 = 'Astrocyte')


# To save your DE results:
write.csv(x = my_de, file = 'my_de.csv', quote = FALSE)



# DE results visualization ------------------------------------------------

# Classic DE test visualization is a volcano plot.

# Bare bones volcano plot
my_volcano <- VolcanoPlot(de_results = my_de)
my_volcano

# Label all genes above a certain log fold-change threshold
my_volcano <- VolcanoPlot(de_results = my_de, label_by_fc = 2)
my_volcano

# Label all genes below a given p-value
my_volcano <- VolcanoPlot(de_results = my_de, label_by_pval = 1e-150)
my_volcano

# Joint strategy to filter by both.
my_volcano <- VolcanoPlot(de_results = my_de, 
                          label_by_fc = 1,
                          label_by_pval = 1e-150,
                          label.size = 5)
my_volcano

# To return the list of genes that meet the above thresholds (pvalue and logFC),
# specify "return_result" to TRUE. NOTE: This will not produce a volcano plot. 
# It will only return the data.frame with the filtered results ie. you will have
# to run VolcanoPlot() twice with "return_result = TRUE" and "return_result = FALSE")
my_volcano <- VolcanoPlot(de_results = my_de,
                          label_by_fc = 1,
                          label_by_pval = 1e-100,
                          return_result = TRUE)
