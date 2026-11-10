load("D:/RProject/2025NBL/sce_all6.rds")
library(Seurat)
library(stringr)
library(harmony)
library(ggplot2)
library(DropletUtils)#Empty-droplet analysis
library(scDblFinder)#Doublet detection
library(scuttle)
library(irlba)
library(sctransform)
library(glmGamPoi)
sce_all6

table(sce_all6$new_tcell_cluster)
DimPlot(sce_all6)
# Create a logical index retaining the target cell populations
keep_cells <- !sce_all6$new_tcell_cluster %in% c('clusterpoliferative_T_cell','cluster16',"cluster19", "cluster20")

# Subset the single-cell object
sce_all6 <- sce_all6[, keep_cells]


sce_all6$newthree=sce_all6$new_tcell_cluster
sce_all6$newthree=gsub('cluster6|cluster11|cluster14','MiddleRisk',sce_all6$newthree)
sce_all6$newthree=gsub('cluster3|cluster4|cluster7|cluster18|cluster16','HighRisk',sce_all6$newthree)
sce_all6$newthree=gsub('cluster10|cluster12|cluster13|cluster15|cluster17','LowRisk',sce_all6$newthree)
sce_all6$newthree=gsub('cluster10|cluster12|cluster13|cluster15|cluster17|cluster0|cluster1|cluster2|cluster8|cluster5|cluster9','LowRisk',sce_all6$newthree)



# Check the resulting annotations
table(sce_all6$newthree)
sce_all6$newthree=gsub('cluster','',sce_all6$newthree)
table(sce_all6$newthree)


library(CellChat)

library(CellChat)
library(Seurat)
library(patchwork)


{
  # Define the cell-type display order
  celltype_levels <- c(
    # Tumor cells first
    "HighRisk","MiddleRisk", "LowRisk", 
    # Stromal and structural cells
    "Endothelial cells", "Epithelial", "Fibroblast",
    # Myeloid cells, from monocytes to differentiated subtypes
    "Mo", "Myeloid cells", "M1", "M2", "Mg_TAM",
    
    # T-cell populations
    "CD8_T_cell", "naive_T_cell", "T cells", "T_regulatory_cell",
    
    # B-cell lineage
    "B cells", "Plasma cells",
    # Other immune cells
    "DC", "NK cells"
  )
  
  # Convert cell-type metadata to a factor with specified levels
  sce_all6@meta.data$newthree <- factor(
    x = sce_all6@meta.data$newthree,
    levels = celltype_levels
  )
  
  # Check factor-level order
  levels(sce_all6$newthree)
  table(sce_all6$newthree)
}
# Inspect cell annotations in the Seurat object
# Check that metadata contains cell-type annotations
print(head(sce_all6@meta.data))

# Extract the expression matrix and cell annotations from SCT or RNA
# This workflow uses SCT-normalized data after SCTransform
data.input <- GetAssayData(sce_all6, assay = "SCT", slot = "data")  # Use SCT-normalized data
# data.input <- GetAssayData([sce_all6], assay = "RNA", slot = "data") # Alternatively use the RNA assay

meta <- sce_all6@meta.data # Extract metadata
meta$newthree <- factor(meta$newthree, levels = celltype_levels)

cell.use <- rownames(meta) # Extract cell identifiers; optionally restrict to a subset

# Create the CellChat object
cellchat <- createCellChat(object = data.input[, cell.use],
                           meta = meta,
                           group.by = "newthree")

# Use the built-in human or mouse CellChat database
CellChatDB <- CellChatDB.human  # Human database
# CellChatDB <- CellChatDB.mouse # Mouse database

# Optionally select a database subset such as Secreted Signaling or Cell-Cell Contact
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use

# Preprocess and identify overexpressed genes
cellchat <- subsetData(cellchat)  # Required preprocessing step
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probabilities; this may take time
cellchat <- computeCommunProb(cellchat, 
                              type = "truncatedMean", # Use a truncated mean
                              trim = 0.1, # Trim the upper and lower 10 percent
                              population.size = TRUE) # Account for population sizes

# Optionally filter low-confidence communication
cellchat <- filterCommunication(cellchat, min.cells = 10) # Require at least 10 cells

# Compute pathway-level communication
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate the communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(meta$celltype))


par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F)

#Heatmap
netVisual_heatmap(cellchat, 
                  measure = "count",  # Alternatively use weight for interaction strength
                  color.heatmap = "Spectral",
                  title.name = "Interaction Counts/Weights")

pathway=cellchat@netP$pathways

dev.off()

length(pathway)

netVisual_chord_gene(cellchat, 
                     signaling = pathway[1],
                     lab.cex = 0.8,
                     legend.pos.y = 30,
                     title.name = paste0(pathway[1], " Ligand-Receptor Pairs"))


# Set the output PDF path
pdf_file <- "D:/RProject/2025NBL/cellchat/All_Pathways_Chord_Plots.pdf"

# Open the PDF device
pdf(pdf_file, width = 12, height = 10)  # Set page dimensions

# Iterate over all pathways
for (sig_pathway in pathway) {
  # Report progress
  cat("\nProcessing pathway:", sig_pathway, "\n")
  
  # Catch errors so one pathway does not interrupt the loop
  tryCatch({
    # Draw the chord diagram
    netVisual_chord_gene(cellchat,
                         signaling = sig_pathway,
                         lab.cex = 0.6,        # Reduce label size
                         small.gap = 1,        # Small gap between adjacent sectors
                         big.gap = 5,          # Gap between groups
                         link.lwd = 0.5,      # Link line width
                         title.name = paste0(sig_pathway, " L-R Pairs"))
    
    # Start a new plot page
    plot.new()
    
  }, error = function(e) {
    # Report the error
    message(paste0("Error in pathway ", sig_pathway, ": ", e$message))
  })
}

# Close the PDF device
dev.off()
{
  
  pdf_file <- "D:/RProject/2025NBL/cellchat/All_Pathways_Chord_Plots1.3.pdf"
  
  # Open the PDF device
  pdf(pdf_file, width = 5, height = 5)  # Set page dimensions
  
  # Iterate over all pathways
  for (sig_pathway in pathway) {
    # Report progress
    cat("\nProcessing pathway:", sig_pathway, "\n")
    
    # Catch errors so one pathway does not interrupt the loop
    tryCatch({
      # Draw the chord diagram
      netVisual_chord_cell(cellchat,
                           signaling = sig_pathway,
                           lab.cex = 0.6,        # Reduce label size
                           small.gap = 1,        # Small gap between adjacent sectors
                           big.gap = 5,          # Gap between groups
                           link.lwd = 0.5,      # Link line width
                           title.name = paste0(sig_pathway, " L-R Pairs"))
      
      # Start a new plot page
      plot.new()
      
    }, error = function(e) {
      # Report the error
      message(paste0("Error in pathway ", sig_pathway, ": ", e$message))
    })
  }
  
  # Close the PDF device
  dev.off()
  
}


netVisual_chord_cell(cellchat,
                     signaling = pathway[8],
                     lab.cex = 0.6,        # Reduce label size
                     small.gap = 1,        # Small gap between adjacent sectors
                     big.gap = 5,          # Gap between groups
                     link.lwd = 0.5,      # Link line width
                     title.name = paste0(pathway[8], " L-R Pairs"))
