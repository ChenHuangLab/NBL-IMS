## 
### ---------------
###
### Create: Jianming Zeng
### Date: 2020-02-09 16:46:35
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2020-02-09   First version
###
### ---------------



rm(list = ls())   
options(stringsAsFactors = F)
library("ChAMP")
library("minfi")
require(GEOquery)
require(Biobase)
load(file = 'step1-output.Rdata')

myLoad  
# Time-consuming step; optionally comment out after the first run
if(F){
  myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
  dim(myNorm) 
  pD=myLoad$pd
  save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')
}
load(file = 'step2-champ_myNorm.Rdata')
# Quality control reduces the original 450K probe set
beta.m=myNorm
group_list=myLoad$pd$Group
dim(beta.m) 
# Three standard matrix quality-control plots; original tutorial by Bioinformatics Skill Tree
if(T){
  
  
  dat=t(beta.m)
  dat[1:4,1:4] 
  library("FactoMineR")#Load these packages for PCA visualization
  library("factoextra")  
  # PCA can be slow for hundreds of thousands of methylation probes
  dat.pca <- PCA(dat , graph = FALSE) 
  fviz_pca_ind(dat.pca,
               geom.ind = "point", # show points only (nbut not "text")
               col.ind = group_list, # color by groups
               # palette = c("#00AFBB", "#E7B800"),
               addEllipses = TRUE, # Concentration ellipses
               legend.title = "Groups"
  )
  ggsave('all_samples_PCA.png')
  
  dat=beta.m
  dat[1:4,1:4] 
  cg=names(tail(sort(apply(dat,1,sd)),1000))#Compute row standard deviations, sort, and select the top 1000 probes
  library(pheatmap)
  pheatmap(dat[cg,],show_colnames =F,show_rownames = F) #Plot the matrix restricted to the selected 1000 probes
  n=t(scale(t(dat[cg,]))) # Scale the selected probe values by row
  n[n>2]=2 
  n[n< -2]= -2
  n[1:4,1:4]
  pheatmap(n,show_colnames =F,show_rownames = F)
  ac=data.frame(group=group_list)
  rownames(ac)=colnames(n)  
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=ac,filename = 'heatmap_top1000_sd.png')
  dev.off()
  
  exprSet=beta.m
  pheatmap::pheatmap(cor(exprSet)) 
  # Inspect whether within-group similarity exceeds between-group similarity
  colD=data.frame(group_list=group_list)
  rownames(colD)=colnames(exprSet)
  pheatmap::pheatmap(cor(exprSet),
                     annotation_col = colD,
                     show_rownames = F,
                     filename = 'cor_all.png')
  dev.off() 
  exprSet=exprSet[names(sort(apply(exprSet, 1,mad),decreasing = T)[1:500]),]
  dim(exprSet)
  # M=cor(log2(exprSet+1)) 
  M=cor(exprSet)
  pheatmap::pheatmap(M,annotation_col = colD)
  pheatmap::pheatmap(M,
                     show_rownames = F,
                     annotation_col = colD,
                     filename = 'cor_top500.png')
  dev.off() 
  
}


myLoad 
beta.m=myLoad$beta
# Legacy alternative normalization workflow using wateRmelon
if(F){
  # Normalize using wateRmelon
  library("wateRmelon")
  beta.m=beta.m[rowMeans(beta.m)>0.005,]
  pdf(file="rawBox.pdf")
  boxplot(beta.m,col = "blue",xaxt = "n",outline = F)
  dev.off()
  beta.m = betaqn(beta.m)
  pdf(file="normalBox.pdf")
  boxplot(beta.m,col = "red",xaxt = "n",outline = F)
  dev.off()
  
  # Then perform basic quality control
  group_list=myLoad$pd$Group
  pdf(file="densityBeanPlot.pdf")
  par(oma=c(2,10,2,2))
  densityBeanPlot(beta.m, sampGroups = group_list)
  dev.off()
  pdf(file="mdsPlot.pdf")
  mdsPlot(beta.m, numPositions = 1000, sampGroups = group_list)
  dev.off()
  
  # Perform downstream differential analysis on beta.m, for example using minfi
  grset=makeGenomicRatioSetFromMatrix(beta.m,what="Beta")
  M = getM(grset)
  # Statistical tests can be slow for hundreds of thousands of methylation probes
  dmp <- dmpFinder(M, pheno=group_list, type="categorical")
  dmpDiff=dmp[(dmp$qval<0.05) & (is.na(dmp$qval)==F),]
  dim(dmpDiff)
}


