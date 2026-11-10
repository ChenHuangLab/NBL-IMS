load("D:/RProject/2025NBL/sce_all6.rds")
sce_all6
cells_to_keep <- sce_all6$newthree != "poliferative_T_cell"
cells_to_keep <- !(sce_all6$new_tcell_cluster %in% c("cluster16",'cluster10','cluster17', "cluster19", "cluster20",'clusterpoliferative_T_cell')) 
# Create the filtered Seurat object
sce_all6 <- subset(sce_all6, cells = colnames(sce_all6)[cells_to_keep])

sce_all6$newthree=sce_all6$new_tcell_cluster
sce_all6$newthree=gsub('cluster6|cluster11|cluster14|cluster18','MiddleRisk',sce_all6$newthree)
sce_all6$newthree=gsub('cluster3|cluster4|cluster7','HighRisk',sce_all6$newthree)
sce_all6$newthree=gsub('cluster0|cluster1|cluster2|cluster8|cluster5|cluster9|cluster12|cluster13|cluster15','LowRisk',sce_all6$newthree)
sce_all6$newthree=gsub('cluster','',sce_all6$newthree)
gc()
table(sce_all6$newthree)

# Create broad cell-type annotations
sce_all6@meta.data <- sce_all6@meta.data %>%
  mutate(celltype_broad = case_when(
    # T-cell class
    newthree %in% c("CD8_T_cell", "naive_T_cell", "T cells", "T_regulatory_cell") ~ "T cells",
    
    # B-cell class
    newthree %in% c("B cells", "Plasma cells") ~ "B lineage",
    
    # Myeloid class
    newthree %in% c("DC", "M1", "M2", "Mg_TAM", "Mo", "Myeloid cells") ~ "Myeloid",
    
    # NK cells
    newthree == "NK cells" ~ "NK cells",
    
    # Stromal cells
    newthree %in% c("Endothelial cells", "Epithelial", "Fibroblast") ~ "Stromal",
    
    # Tumor-cell risk states
    newthree %in% c("HighRisk", "LowRisk", "MiddleRisk") ~ "Tumor cells",
    
    # Label unclassified cells as Other
    TRUE ~ "Other"
  ))

# Check broad cell-type annotations
table(sce_all6$celltype_broad)

# Example output:
# B lineage      Myeloid      NK cells      Stromal    T cells  Tumor cells 
#        4244        14695         7407        12720        19667       118779 

# Optionally extract broad cell-type subsets
# Extract all T cells
t_cells <- subset(sce_all6, subset = celltype_broad == "T cells")

# Extract all myeloid cells
myeloid <- subset(sce_all6, subset = celltype_broad == "Myeloid")

# Extract macrophage subpopulations for finer annotation
macrophages <- subset(sce_all6, 
                      subset = newthree %in% c("M1", "M2", "Mg_TAM"))
DimPlot(macrophages,group.by = 'newthree')
DimPlot(t_cells,group.by = 'newthree')

# Normalization and PCA
macrophages <- NormalizeData(macrophages) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(npcs = 30)

# Batch adjustment using Harmony
macrophages <- RunHarmony(
  macrophages,
  group.by.vars = "metaset",  # Batch metadata column
  assay.use = "RNA",
  plot_convergence = TRUE
)

# Dimension reduction and clustering
macrophages <- RunUMAP(macrophages, reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)
DimPlot(macrophages,group.by = 'newthree')
Idents(macrophages)=macrophages$newthree
macrophages.marker=FindAllMarkers(macrophages,all.pos=T)


# Normalization and PCA


# Batch adjustment using Harmony
t_cells <- RunHarmony(
  t_cells,
  group.by.vars = "metaset",
  assay.use = "RNA",
  plot_convergence = TRUE
)

# Dimension reduction and clustering
t_cells <- RunUMAP(t_cells, reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)
Idents(t_cells)=t_cells$newthree
t_cells.marker=FindAllMarkers(t_cells,all.pos=T)
DimPlot(t_cells,group.by = 'newthree')

# Load required packages
library(dplyr)
library(Seurat)
library(ggplot2)

# Select the top three marker genes per cluster by avg_log2FC
top3_genes <- t_cells.marker %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene)
# Original marker note: naive CCR7,TCH7,LEF1
# Create a vertical dot plot
# Specify a custom display order
custom_order <- c("CD8_T_cell", "T cells", "T_regulatory_cell", "naive_T_cell")
t_cells$newthree <- factor(t_cells$newthree, levels = custom_order)

DotPlot(t_cells, 
        features = top3_genes,
        group.by = "newthree",
        cols = c("lightgrey", "red"),
        dot.scale = 8) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) # Rotate x-axis labels

# Display the plot
print(dotplot)
unique(macrophages$newthree)
custom_order1 <- c("Mg_TAM", "M1", "M2")
macrophages$newthree <- factor(macrophages$newthree, levels = custom_order1)

unique(macrophages$newthree)
top3_genes1 <- macrophages.marker %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC)) %>%
  pull(gene)
DotPlot(macrophages, 
        features = top3_genes1,
        group.by = "newthree",
        cols = c("lightgrey", "red"),
        dot.scale = 8) + 
  coord_flip() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) # Rotate x-axis labels
save(macrophages.marker,t_cells.marker,file = 't_cells.mo.marker.rds')
