load('/hdd/1/fengyunwen/rds/ne2.1.rds')
library(Seurat)


Idents(ne2_1)=ne2_1$newthree

name=FindAllMarkers(ne2_1,only.pos = T)
name1=name$gene[which(name$avg_log2FC>=1)]

Idents(ne2_1)=ne2_1$three
orig3.marker=FindAllMarkers(ne2,only.pos = T)



load('/hdd/1/fengyunwen/rds/ne2.1.rds')
ne2_1$newthree1=ne2_1$newthree
ne2_1$newthree1=gsub('good','HighRisk',ne2_1$newthree1)
ne2_1$newthree1=gsub('bad','LowRisk',ne2_1$newthree1)
ne2_1$newthree1=gsub('mixed','MiddleRisk',ne2_1$newthree1)
Idents(ne2_1)=ne2_1$newthree1
newthree1.marker=FindAllMarkers(ne2_1,only.pos = T)
save(newthree1.marker,file = '/hdd/1/fengyunwen/rds/newthree1.marker.rds')


Idents(ne2)=ne2$three
orig3.marker=FindAllMarkers(ne2,only.pos = T)



type_exp <- AverageExpression(
  ne2_1[name1,],
  assays = "RNA",      # Use the RNA assay
  slot = "counts",       # Use the selected expression slot; counts contains raw counts
  group.by = "newthree"
)$RNA

data_exp <- AverageExpression(
  ne2_1,
  assays = "RNA",     # Use the RNA assay
  slot = "data",      # Use the selected normalized expression slot
  group.by = "newthree"
)$RNA

data_exp1 <- AverageExpression(
  ne2_1,
  assays = "SCT",     # Use the RNA assay
  slot = "scale.data",      # Use the selected normalized expression slot
  group.by = "newthree"
)$RNA


data_exp11 <- AverageExpression(
  ne2_1,
  assays = "SCT",     # Use the RNA assay
  slot = "data",      # Use the selected normalized expression slot
  group.by = "newthree"
)$SCT

counts.exp <- AverageExpression(
  ne2_1,
  assays = "RNA",      # Use the RNA assay
  slot = "counts",       # Use the selected expression slot; counts contains raw counts
  group.by = "newthree"
)$RNA

max(AverageExpression(
  ne2_1,
  assays = "RNA",      # Use the RNA assay
  slot = "counts",       # Use the selected expression slot; counts contains raw counts
  group.by = "newthree"
)$RNA)

inf.exp=AverageExpression(
  ne2_1,
  assays = "RNA",      # Use the RNA assay
  slot = "data",       # Use the selected expression slot; counts contains raw counts
  group.by = "newthree"
)$RNA


save(counts.exp,file = '/hdd/1/fengyunwen/rds/counts.exp.rds')

save(inf.exp,file = '/hdd/1/fengyunwen/rds/inf.exp.rds')

save(type_exp,file = '/hdd/1/fengyunwen/rds/newthree_exp.rds')

save(data_exp11,file = '/hdd/1/fengyunwen/rds/data_exp11.rds')

load("D:/RProject/2025NBL/generated_data/data_exp11.rds")


write.table(type_exp,file = '/hdd/1/fengyunwen/csv/newthree_exp.txt')



load("D:/RProject/2025NBL/generated_data/bulk_and_external.rds")
load("D:/RProject/2025NBL/generated_data/newthree_exp.rds")

library(purrr)      # Provides purrr::reduce
library(IRanges)    # Provides IRanges::reduce
library(GenomicRanges) # Provides GenomicRanges::reduce
library(e1071)
library(preprocessCore)
library(parallel)
library(bseqsc)#Includes many CIBERSORT dependencies; may help resolve installation dependencies
library(CIBERSORT)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)

load('D:/RProject/2025NBL/generated_data/group1.rds')


results.with16 <- cibersort(
  sig_matrix = type_exp,
  mixture_file = target_new,
  perm = 1000  # Number of permutations
)

data_exp11

results.sct <- cibersort(
  sig_matrix = data_exp11,
  mixture_file = target_new,
  perm = 1000  # Number of permutations
)


ne3=data.frame(results.sct)
results.inf.new3

ne3$sample <- rownames(ne3)  # Add sample identifiers
ne3 <- merge(ne3, group1[,c('X','group')], by.x = "sample", by.y = "X")  # Merge group annotations
ne3=ne3[!duplicated(ne3$sample),]


colMeans(ne3[which(ne3$group%in%'NMF_1'),c(2,3,4)])
colMeans(ne3[which(ne3$group%in%'NMF_2'),c(2,3,4)])
colMeans(ne3[which(ne3$group%in%'NMF_3'),c(2,3,4)])
