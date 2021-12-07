

##### Differential Gene Expression (DE) tests using Seurat #####

# Libraries
library('Seurat')
library('ggplot2')
library('dplyr')
install.packages('ggrepel') # if you do not have 'ggrepel' installed.
library('ggrepel')

# Set working directory
setwd(dir = 'D:/introductory_scSeq/')

# Load data
macroglia <- readRDS(file = './data/macroglia.rds')

# DE tests _must_ be done with the RNA assay. 
DefaultAssay(macroglia) <- 'RNA'

# DE tests compare between identities
Idents(macroglia) <- 'macroglia_7dpi'



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


# You can also split a certain cell-type by time and perform DE between time 
# points.
#   1) Specify the cell-type(s) with "subset.ident"
#   2) Set "group.by" to "time"
#   3) Select the two groups (time-points)
my_de <- FindMarkers(macroglia,
                     subset.ident = 'Astrocyte',
                     group.by = 'time',
                     ident.1 = '3dpi',
                     ident.2 = '1dpi')
my_de <- FindMarkers(macroglia,
                     subset.ident = c('Ependymal-A','Ependymal-B','Ependymal-C'),
                     group.by = 'time',
                     ident.1 = '7dpi',
                     ident.2 = '3dpi')


# Conversely, take all cells from a time point, then group cells by subcluster
# identities, finally compare two groups.
Idents(macroglia) <- 'time'
my_de <- FindMarkers(macroglia,
                     subset.ident = '1dpi',
                     group.by = 'macroglia_subcluster',
                     ident.1 = 'Astroependymal',
                     ident.2 = 'Astrocyte')


# To save your DE results:
write.csv(x = my_de, file = './results/my_de.csv', quote = FALSE)

# Re-load DE results
my_de <- read.csv(file = './results/my_de.csv', row.names = 1)



# DE results visualization ------------------------------------------------

# Classic DE test visualization is a volcano plot.

# Run the following to create a function that will generate a volcano plot. It
# will take as input the data.frame that is generated from the FindMarkers() 
# function.
VolcanoPlot <- function(
  de_results, 
  label_by_fc = NULL, 
  label_by_pval = NULL,
  label_size = 4
) {
  require('ggrepel')
  de_results <- de_results %>%
    mutate(
      'gene' = rownames(.),
      'log_pval' = -log10(p_val_adj),
      'label' = ifelse(test = abs(avg_log2FC) <= label_by_fc &
                         p_val_adj >= label_by_pval,
                       yes = NA, 
                       no = gene)
    )
  # return(de_results)
  de_volcano <- de_results %>%
    ggplot(mapping = aes(x = avg_log2FC, y = log_pval)) +
    geom_point(color = 'grey40') + 
    geom_vline(xintercept = 0, linetype = 'dashed', size = 1, color = 'black') +
    scale_x_continuous() +
    xlab(label = 'Log2(fold-change)') +
    ylab(label = '-Log10(adjusted p-value)') +
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          panel.border = element_rect(fill = NA, color = 'black', size = 1),
          axis.text = element_text(size = 14, color = 'black'),
          axis.title = element_text(size = 16, color = 'black'),
          panel.grid = element_line(color = 'grey60', linetype = 'dotted')) +
    geom_text_repel(mapping = aes(label = label), size = label_size) +
    geom_vline(xintercept = label_by_fc, color = 'red', linetype = 'dashed',
               size = 1) +
    geom_vline(xintercept = -label_by_fc, color = 'red', linetype = 'dashed',
               size = 1) +
    geom_hline(yintercept = -log10(label_by_pval), color = 'blue', 
               linetype = 'dashed', size = 1)
  return(de_volcano)
}



# Volcano plot
my_volcano <- VolcanoPlot(de_results = my_de, 
                          label_by_fc = 1,
                          label_by_pval = 1e-50)
my_volcano
