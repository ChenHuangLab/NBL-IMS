setwd('/home/yuwenfeng/NBLdata/7.5/')

{library(Seurat)
  library(stringr)
  library(harmony)
  library(ggplot2)
  library(DropletUtils)#Empty-droplet analysis
  library(scDblFinder)#Doublet detection
  library(scuttle)
  library(irlba)
  library(sctransform)
  library(glmGamPoi)
}
dir_path='/home/yuwenfeng/NBLdata/GSE140819/h5/h5/'
dir_file=list.files(dir_path)
file_list=paste0(dir_path,dir_file)

nb1=list()
nb1 <-  lapply(file_list, function(x){
  
  a=Read10X_h5( x )
  a[1:4,1:4]
  sce <- CreateSeuratObject(counts = a,project = gsub(".*GSM([0-9]+).*", "GSM\\1", x), min.cells = 3, min.features = 200 )
})

sce_GSE140819=merge(nb1[[1]],nb1[-1])


dir_path='/home/yuwenfeng/NBLdata/GSE192906/'
dir_file=list.files(dir_path)
file_list=paste0(dir_path,dir_file)

nb1=list()
nb1 <-  lapply(file_list, function(x){
  
  a=read.delim( x )
  a[1:4,1:4]
  sce <- CreateSeuratObject(counts = a,project = gsub(".*GSM([0-9]+).*", "GSM\\1", x), min.cells = 3, min.features = 200 )
})
sce1=list()
nb1[[1]]$dataset='GSM5768743'
nb1[[2]]$dataset='GSM5768744'
nb1[[3]]$dataset='GSM5768745'
nb1[[4]]$dataset='GSM5768746'
nb1[[5]]$dataset='GSM5768747'


sce_GSE192906=merge(nb1[[1]],nb1[-1])

sce_GSE140819$dataset=sce_GSE140819$orig.ident


table(sce_GSE192906$dataset)
table(sce_GSE140819$dataset)


save(sce_GSE140819,sce_GSE192906,file = '/home/yuwenfeng/NBLdata/sce_GSE140819__sce_GSE192906-new.rds')

load('/home/yuwenfeng/NBLdata/db_all.rds')
load('/home/yuwenfeng/NBLdata/sce_GSE140819__sce_GSE192906-new.rds')
nb_GOSH=readRDS('/home/yuwenfeng/NBLdata/nb_GOSH.rds')

sce_all = list()
db_all$dataset=db_all$orig.ident

sce_all[["GSE137804"]] = db_all 
sce_all[["GSE140819"]] = sce_GSE140819
sce_all[["GSE192906"]] = sce_GSE192906

nb_GOSH=UpdateSeuratObject(nb_GOSH)
sce_all[["nb_GOSH"]]= nb_GOSH

sce_all$GSE137804$metaset='GSE137804'
sce_all$GSE140819$metaset='GSE140819'
sce_all$GSE192906$metaset='GSE192906'
sce_all$nb_GOSH$metaset='nb_GOSH'
sce_all$nb_GOSH$dataset='nb_GOSH'

print("All datasets have been list !")


###Remove doublets
for(metaset in names(sce_all)){
  print(metaset)
  samples <- unique(sce_all[[metaset]]$dataset)
  ####Remove doublets
  ###Identify doublets
  doublet_scores <- c()
  doublet_predictions <- c()
  for(sample in samples) {
    print(sample)
    # Extract sample data
    sce <- sce_all[[metaset]][,sce_all[[metaset]]$dataset == sample]
    # Extract the expression matrix
    counts <- sce@assays$RNA@counts
    # Run scDblFinder
    res <- scDblFinder(counts)
    # Add results to the Seurat object
    doublet_scores_tmp <- res$scDblFinder.score
    doublet_predictions_tmp <- res$scDblFinder.class
    # Store results in the Seurat object
    doublet_scores <- c(doublet_scores,doublet_scores_tmp)
    doublet_predictions <- c(doublet_predictions,as.character(doublet_predictions_tmp))
  }
  sce_all[[metaset]][["doublet_scores"]] <- doublet_scores
  sce_all[[metaset]][["doublet_predictions"]] <- doublet_predictions
}

for (metaset in names(sce_all)) {
  print(metaset)
  sce_all[[metaset]] <- subset(sce_all[[metaset]], subset = doublet_predictions != "doublet")
}



###Quality control
for (metaset in names(sce_all)) {
  print(metaset)
  sce_all[[metaset]][["percent.mt"]] <- PercentageFeatureSet(sce_all[[metaset]], pattern = "^MT-")
  HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
  HB_m <- match(HB.genes, rownames(sce_all[[metaset]]@assays$RNA)) 
  HB.genes <- rownames(sce_all[[metaset]]@assays$RNA)[HB_m] 
  HB.genes <- HB.genes[!is.na(HB.genes)] 
  sce_all[[metaset]][["percent.HB"]]<-PercentageFeatureSet(sce_all[[metaset]], features=HB.genes) 
  
  samples <- unique(sce_all[[metaset]]$dataset)
  nFeature_RNA_filter <- c()
  nCount_RNA_filter <- c()
  for(sample in samples){
    print(sample)
    dat <- sce_all[[metaset]]@meta.data[which(sce_all[[metaset]]@meta.data$dataset == sample ),]
    nFeature_RNA_filter_tmp <- isOutlier(dat$nFeature_RNA, nmads=4, type="both", log=FALSE)
    nCount_RNA_filter_tmp <- isOutlier(dat$nCount_RNA, nmads=4, type="both", log=FALSE)
    nFeature_RNA_filter <- c(nFeature_RNA_filter,nFeature_RNA_filter_tmp)
    nCount_RNA_filter <- c(nCount_RNA_filter,nCount_RNA_filter_tmp)
  }
  sce_all[[metaset]][["nFeature_RNA_filter"]] <- nFeature_RNA_filter
  sce_all[[metaset]][["nCount_RNA_filter"]] <- nCount_RNA_filter
}


sce_all1 <- list()
for (metaset in names(sce_all)) {
  sce_all1[[metaset]] <- subset(sce_all[[metaset]], 
                                subset = nFeature_RNA_filter == FALSE & 
                                  nCount_RNA_filter == FALSE &
                                  percent.mt < 10 &
                                  percent.HB < 3)
}

print("QC have done !")
save(sce_all1,file = "sce_all1.RData")


###SCTransform-based integration and batch adjustment
sce_all2 <- lapply(sce_all1,function(object) {
  object <- SCTransform(object,verbose = FALSE,vst.flavor="v2",method = "glmGamPoi",vars.to.regress = "percent.mt")
})

integration_features <- SelectIntegrationFeatures(object.list = sce_all2,nfeatures = 3000)
sce_all2 <- PrepSCTIntegration(object.list = sce_all2, anchor.features = integration_features)
integration_anchors <- FindIntegrationAnchors(object.list = sce_all2, normalization.method = "SCT", anchor.features = integration_features)

sce_all2 <- IntegrateData(anchorset = integration_anchors,normalization.method = "SCT")
#Use the integrated assay as the default for downstream analysis
DefaultAssay(sce_all2) <- "integrated"
sce_all2 <- RunPCA(object = sce_all2, npcs = 30, verbose = FALSE)
sce_all2 <- FindNeighbors(object = sce_all2, reduction = "pca", dims = 1:30)
sce_all2 <- RunUMAP(object = sce_all2, reduction = "pca", dims = 1:30)
sce_all2 <- RunTSNE(sce_all2, dims = 1:30, reduction = "pca",check_duplicates = FALSE)

print("scTransform Done!")
save(sce_all2,file = "sce_all2.RData")


###
pdf("UMAP_plot_dataset_sctf.pdf",height = 8, width = 8)
DimPlot(sce_all2, 
        reduction = "umap",
        group.by = "dataset", 
        pt.size=0.5, 
        label = F,
        repel = TRUE, 
        raster = FALSE)
dev.off()


pdf("UMAP_plot_metadata_sctf.pdf",height = 8, width = 8)
DimPlot(sce_all2, 
        reduction = "umap",
        group.by = "metaset", 
        pt.size=0.5, 
        label = F,
        repel = TRUE, 
        raster = FALSE)
dev.off()

save(sce_all2,file = "sce_all2.RData")

###Inspect clustering resolutions
library(clustree)
obj <- FindClusters(sce_all2, resolution = seq(0.1,1,by=0.1))
pdf("clustree.pdf",height =10, width = 10)
clustree(obj)
dev.off()

sce_all2 <- FindClusters(sce_all2, resolution = 0.5)
sce_all2 <- RunUMAP(object = sce_all2, reduction = "pca", dims = 1:30)
sce_all2 <- RunTSNE(sce_all2, dims = 1:30, reduction = "pca",check_duplicates = FALSE)

pdf("UMAP_plot_umap_sctf.pdf",height = 8, width = 8)
DimPlot(sce_all2, 
        reduction = "umap", 
        pt.size=0.5, 
        label = T,
        repel = TRUE,
        raster = FALSE)
dev.off()

pdf("UMAP_plot_tsne_sctf.pdf",height = 8, width = 8)
DimPlot(sce_all2, 
        reduction = "tsne", 
        pt.size=0.5, 
        label = T,
        repel = TRUE,
        raster = FALSE)
dev.off()



a1 <- DimPlot(sce_all2, reduction = "umap",group.by = "seurat_clusters", pt.size=0.5, label = F,repel = TRUE, raster = FALSE)
a2 <- DimPlot(sce_all2, reduction = "umap",group.by = "metaset", pt.size=0.5, label = F,repel = TRUE, raster = FALSE)
a3 <- DimPlot(sce_all2, reduction = "umap",group.by = "dataset", pt.size=0.5, label = F,repel = TRUE, raster = FALSE)

pdf("seuratclusters_dataset_origident.pdf",onefile = F,height = 10,width = 30)
a1+a2+a3
dev.off()

# Find marker genes across all clusters
sce_all2 <- PrepSCTFindMarkers(sce_all2, assay = "SCT")

cluster_markers <- FindAllMarkers(object = sce_all2, 
                                  assay = "SCT",
                                  only.pos = TRUE, 
                                  min.pct = 0.25, 
                                  logfc.threshold = 0.25)
write.csv(cluster_markers,"cluster_markers.csv")
integrated_snn_res.0.5

###Cell-type marker genes
marker_genes <- list('Ne cells'  =c('DBH','SYP','TH','MYCN','NCAM1','PHOX2A','PHOX2B','NYAP2'),
                     'T cells' = c("CD3D","CD8A","IL7R"),#CD3E,CD2,CD3G
                     'NK cells' = c("KLRD1","KLRF1","GNLY","CX3CR1","XCL2"),# FGFBP2 has low expression; FCG3RA is not expressed
                     'B cells' = c("CD79A","CD19","MS4A1"),#IGHG1
                     'Myeloid cells' = c("CXCL8","S100A8","S100A9"),
                     'Fibroblast' = c("COL1A1","COL1A2","DCN","FGF7"),
                     'Epithelial' = c("EPCAM","KRT19","ALDH1A1", "CD24"),                  
                     'Macrophage' = c("CD68","MSR1","APOE","TREM2","FOLR2"),
                     'Endothelial cells' = c("PECAM1","CLDN5","VWF"),
                     'Plasma cells' = c("MZB1"))##Cell marker heatmap


##Cell marker dot plot
pdf("dotplot_marker_genes_new.pdf",onefile = F,height = 8,width = 12)
DotPlot(sce_all2,
        assay = "SCT",
        features = marker_genes,group.by = 'seurat_clusters')+
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=8), #Set x-axis tick label angle, justification, font, and size
        axis.text.y=element_text(family="Times",size=8,face="plain"), #Set y-axis tick label font, size, and style
        axis.title.x=element_blank(), #Set x-axis title properties
        axis.title.y=element_blank(), #Set y-axis title properties
        legend.text=element_text(face="italic", family="Times", colour="black",  #Set legend label text properties
                                 size=8),
        legend.title=element_text(face="italic", family="Times", colour="black", #Set legend title properties
                                  size=8),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #Remove the panel border and style the axis lines
        panel.grid.major = element_blank(),   #Hide grid lines
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=7))  #Hide grid lines

dev.off()

pdf("dimplot_marker_seurat_clusters.pdf",onefile = F,height = 8,width = 8)
DimPlot(sce_all2,
        reduction = "umap",group.by = 'seurat_clusters',
        pt.size=0.5, 
        label = F,
        repel = TRUE, 
        raster = FALSE,cols = c('26'='red'))
dev.off()

pdf("vinnPlot_marker_genes_new.pdf",onefile = F,height = 42,width = 56)
VlnPlot(sce_all2,
        assay = "SCT",
        features = unlist(marker_genes),
        raster = FALSE,
        group.by = "RNA_snn_res.0.4",
        pt.size = 0,ncol = 6)+
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=9), #Set x-axis tick label angle, justification, font, and size
        axis.text.y=element_text(family="Times",size=9,face="plain"), #Set y-axis tick label font, size, and style
        axis.title.x=element_blank(), #Set x-axis title properties
        axis.title.y=element_blank(), #Set y-axis title properties
        legend.text=element_text(face="italic", family="Times", colour="black",  #Set legend label text properties
                                 size=10),
        legend.title=element_text(face="italic", family="Times", colour="black", #Set legend title properties
                                  size=10),
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #Remove the panel border and style the axis lines
        panel.grid.major = element_blank(),   #Hide grid lines
        panel.grid.minor = element_blank())  #Hide grid lines
dev.off()
# Proliferative cells
# ###Manually annotate cluster cell types based on the dot plot
# Idents(sce_all2) <- "seurat_clusters"
# sce_all2 <- RenameIdents(sce_all2,
#                          '0' = 'NE cells',
#                          '1' = 'NE cells',
#                          '3' = 'NE cells',
#                          '6' = 'NE cells',
#                          '7' = 'NE cells',
#                          '22' =  'NE+Epithelial',
#                          "8"= "NK cells",
#                          "9"= "T cells",
#                          '10' = 'NE cells',
#                          '15' = 'NE cells',
#                          '16' = 'NE cells',
#                          "17" = "T cells",
#                          '11' = "T cells",
#                          '18' = "T cells",
#                          "20"= "B cells",
#                          '24' = 'Myeloid cells',
#                          '25' = "T cells",
#                          '37' = 'NK cells',
#                          '38' = 'Fibroblast',
#                          '30' = 'Fibroblast',
#                          '32' = 'NE cells',
#                          '26' = 'NE cells',
#                          '23' = 'NE cells',
#                          '13' = 'Fibroblast',
#                          '12' = 'Macrophage',
#                          '19' = 'Endothelial cells',
#                          '27' = 'NE cells',
#                          '29' = 'Macrophage',
#                          '31' = 'NE cells',
#                          '34' = 'Plasma cells',
#                          '35' = 'Epithelial',
#                          '33' = 'Epithelial',
#                          '36' = 'Plasma cells',
#                          '25' =  'NE+Epithelial',
#                          '14' =  'NE+Epithelial',
#                          '5' =  'NE+Epithelial',
#                          '28' =  'NE+Epithelial',
#                          '2' =  'NE+Epithelial',
#                          '39' = 'special',
#                          '4' = 'special',
#                          '21' = 'NE cells',
#                          '40' = 'NE+Epithelial')
# 
# sce_all2@meta.data$cell_type <- Idents(sce_all2)
# 
# pdf("dotplot_marker_genes_new.celltype.pdf",onefile = F,height = 8,width = 25)
# DotPlot(sce_all2,
#         assay = "SCT",
#         features = marker_genes)+
#   theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=8), #Set x-axis tick label angle, justification, font, and size
#         axis.text.y=element_text(family="Times",size=8,face="plain"), #Set y-axis tick label font, size, and style
#         axis.title.x=element_blank(), #Set x-axis title properties
#         axis.title.y=element_blank(), #Set y-axis title properties
#         legend.text=element_text(face="italic", family="Times", colour="black",  #Set legend label text properties
#                                  size=8),
#         legend.title=element_text(face="italic", family="Times", colour="black", #Set legend title properties
#                                   size=8),
#         panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #Remove the panel border and style the axis lines
#         panel.grid.major = element_blank(),   #Hide grid lines
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_text(size=7))  #Hide grid lines
# dev.off()
# 
# 
# pdf("dimplot_marker_genes_new.celltype.pdf",onefile = F,height = 8,width = 8)
# DimPlot(sce_all2, 
#         reduction = "umap",
#         pt.size=0.5, 
#         label = F,
#         repel = TRUE, 
#         raster = FALSE)
# dev.off()
# 
# 
# group4=c('SST',#Somatostatin
#          'EGR1',
#          'SKA2',
#          'RPS17',
#          'GNB2L1',
#          'FOS',
#          'TTC3')
# 
# Idents(sce_all2) <- "seurat_clusters"
# 
# pdf("dotplot_marker_group4.pdf",onefile = F,height = 8,width = 25)
# DotPlot(sce_all2,
#         assay = "SCT",
#         features = group4)+
#   theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times",size=8), #Set x-axis tick label angle, justification, font, and size
#         axis.text.y=element_text(family="Times",size=8,face="plain"), #Set y-axis tick label font, size, and style
#         axis.title.x=element_blank(), #Set x-axis title properties
#         axis.title.y=element_blank(), #Set y-axis title properties
#         legend.text=element_text(face="italic", family="Times", colour="black",  #Set legend label text properties
#                                  size=8),
#         legend.title=element_text(face="italic", family="Times", colour="black", #Set legend title properties
#                                   size=8),
#         panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), #Remove the panel border and style the axis lines
#         panel.grid.major = element_blank(),   #Hide grid lines
#         panel.grid.minor = element_blank(),
#         strip.text.x = element_text(size=7))  #Hide grid lines
# dev.off()
# 
# pdf("dimplot_group4.pdf",onefile = F,height = 8,width = 8)
# DimPlot(sce_all2)
# dev.off()
# 
# celltype_markers <- FindAllMarkers(object = sce_all2, 
#                                   only.pos = TRUE, 
#                                   min.pct = 0.25, 
#                                   logfc.threshold = 0.25)
# write.csv(celltype_markers,"celltype_markers.csv")
# sce_all3=sce_all2
# save(sce_all3,'sce_all3.rds')
# 
# 
# ne2=subset(sce_all3,cell_type=c('NE+Epithelial','NE cells'))
# 
# ne2 <- RunPCA(object = ne2, npcs = 30, verbose = FALSE)
# ne2 <- FindNeighbors(object = ne2, reduction = "pca", dims = 1:30)
# ne2 <- RunUMAP(object = ne2, reduction = "pca", dims = 1:30)
# ne2 <- RunTSNE(ne2, dims = 1:30, reduction = "pca",check_duplicates = FALSE)
# 
# 
# ###
# pdf("ne2UMAP_plot_dataset_sctf.pdf",height = 8, width = 8)
# DimPlot(ne2, 
#         reduction = "umap",
#         group.by = "dataset", 
#         pt.size=0.5, 
#         label = F,
#         repel = TRUE, 
#         raster = FALSE)
# dev.off()
# 
# 
# pdf("ne2UMAP_plot_metadata_sctf.pdf",height = 8, width = 8)
# DimPlot(ne2, 
#         reduction = "umap",
#         group.by = "metaset", 
#         pt.size=0.5, 
#         label = F,
#         repel = TRUE, 
#         raster = FALSE)
# dev.off()
# 
# save(ne2,file = "ne2.RData")
# 
# ###Inspect clustering resolutions
# library(clustree)
# objne2 <- FindClusters(ne2, resolution = seq(0.1,1,by=0.1))
# pdf("ne2clustree.pdf",height =10, width = 10)
# clustree(objne2)
# dev.off()
# 
# 
# ne2neepi=subset(sce_all3,cell_type==c('NE+Epithelial'))
# 
# 
# ne2neepi <- RunPCA(object = ne2neepi, npcs = 30, verbose = FALSE)
# ne2neepi <- FindNeighbors(object = ne2neepi, reduction = "pca", dims = 1:30)
# ne2neepi <- RunUMAP(object = ne2neepi, reduction = "pca", dims = 1:30)
# ne2neepi <- RunTSNE(ne2neepi, dims = 1:30, reduction = "pca",check_duplicates = FALSE)
# 
# library(clustree)
# objne2neepi <- FindClusters(ne2neepi, resolution = seq(0.1,1,by=0.1))
# pdf("ne2neepiclustree.pdf",height =10, width = 10)
# clustree(objne2neepi)
# dev.off()
# 
# 
# ne2ne=subset(sce_all3,cell_type==c('NE cells'))
# 
# ne2ne <- RunPCA(object = ne2ne, npcs = 30, verbose = FALSE)
# ne2ne <- FindNeighbors(object = ne2ne, reduction = "pca", dims = 1:30)
# ne2ne <- RunUMAP(object = ne2ne, reduction = "pca", dims = 1:30)
# ne2ne <- RunTSNE(ne2ne, dims = 1:30, reduction = "pca",check_duplicates = FALSE)
# 
# library(clustree)
# objne2ne <- FindClusters(ne2ne, resolution = seq(0.1,1,by=0.1))
# pdf("ne2neclustree.pdf",height =10, width = 10)
# clustree(objne2ne)
# dev.off()
# 
# 
# ne2neepi <- FindClusters(ne2neepi, resolution = 0.1)
# ne2neepi <- RunUMAP(object = ne2neepi, reduction = "pca", dims = 1:30)
# ne2neepi <- RunTSNE(ne2neepi, dims = 1:30, reduction = "pca",check_duplicates = FALSE)
# 
# 
# 
# ne2ne <- FindClusters(ne2ne, resolution = 0.1)
# ne2ne <- RunUMAP(object = ne2ne, reduction = "pca", dims = 1:30)
# ne2ne <- RunTSNE(ne2ne, dims = 1:30, reduction = "pca",check_duplicates = FALSE)
# 
# ne2neepi_marker=FindAllMarkers(ne2neepi, 
#                                only.pos = TRUE, 
#                                min.pct = 0.25, 
#                                logfc.threshold = 0.25)
# 
# 
# 
# ne2ne_marker=FindAllMarkers(ne2ne, 
#                      only.pos = TRUE, 
#                      min.pct = 0.25, 
#                      logfc.threshold = 0.25)
# save(ne2neepi_marker,ne2ne_marker,file='7.10.ne2markers.rds')
# 
# write.csv(ne2neepi_marker,"ne2neepi_marker.csv")
# write.csv(ne2ne_marker,"ne2ne_marker.csv")
# 
# 
# save(ne2ne,file = "ne2ne.RData")
# save(ne2neepi,file = "ne2neepi.RData")
# 
# save(ne2neepi_marker,ne2ne_marker,file='ne2markers.rds')
# 
# rm(list = ls())
# gc()
# load('/home/yuwenfeng/NBLdata/7.10/nb_logrank_test.RData')
# load('ne2markers.rds')
