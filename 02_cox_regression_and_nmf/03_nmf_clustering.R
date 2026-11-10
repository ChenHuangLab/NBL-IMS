rm(list = ls())
options(stringsAsFactors = F)
setwd("D:/RProject/NMF")
load("D:/RProject/metha/10.18cox.RData")
index=intersect(topgene1$gene,colnames(h))
dat=h[,index]
dat=t(dat)
library(NMF)
library(doMPI) #Original workflow uses this package for larger NMF ranks
#install.packages("doMPI")
#Select the NMF rank
# Compare NMF ranks 2 through 7
res <- nmf(dat,2:7,nrun=50)
plot(res)
#nmf(x, rank,method, seed = nmf.getOption("default.seed"),rng = NULL, nrun = if (length(rank) > 1) 30 else 1,model = NULL, .options = list(),pbackend = nmf.getOption("pbackend"),callback = NULL, ...)
#"D:/RProject/NMF"
#Repeat NMF; see the actual rank argument below
seed = 2020820
nmf.rank4 <- nmf(dat[1:16500,], 
                 rank = 2, 
                 nrun=50,
                 seed = 2020820, 
                 method = "brunet")
save=nmf.rank4
ing=dat[1:1000,1:20]
test=nmf(ing[1:1000,], 
         rank = 2, 
         nrun=50,
         seed = 2020820, 
         method = "brunet")
nmf.rank4=test
jco <- c("#2874C5","#EABF00","#C6524A","#868686")
index <- extractFeatures(nmf.rank4,"max") 
sig.order <- unlist(index)

test1=nmf(ing,3,nrun=3)

NMF.Exp.rank4 <- ing[sig.order,]
NMF.Exp.rank4 <- na.omit(NMF.Exp.rank4) #Remove missing rows after feature ordering
group <- predict(nmf.rank4) # Extract predicted subtypes
table(group)
consensusmap(nmf.rank4,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(NMF.Exp.rank4)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))

############################################
library(NMF)
library(Rmpi)

setwd("E:/APA/APA in breast cancer/NMF_clusering_APA")

exp_matrix<-read.table(gzfile("E:/APA/APA in breast cancer/TCGA-dataset/TCGA-BRCA.htseq_fpkm.tsv.gz"),header=T)

rownames(exp_matrix)<-exp_matrix[,1]
exp_matrix<-exp_matrix[,-1]

for(i in 1:length(rownames(exp_matrix))){
  rownames(exp_matrix)[i]<-gsub ("\\.\\d+","",rownames(exp_matrix)[i])
}


APA<-read.table("APA_ensmbel2symbol.uniq.list",header=T)
rownames(APA)=APA[,1]

common<-intersect(APA$Ensembl,rownames(exp_matrix))

exp2<-exp_matrix[common,]

tumor_exp<-exp2[,grep(pattern="01A",names(exp2))]

#Identifier conversion###
cc<-merge(tumor_exp,APA,by="row.names")
rownames(cc)<-cc[,1075]
cc<-cc[,-1]
cc<-cc[,-1074]
cc<-cc[,-1073]
tumor_exp<-cc
write.table(tumor_exp,"APA_exp.txt",quote=F)

res <- nmf(tumor_exp,2:7,nrun=10)

pdf("res.pdf")
plot(res)
dev.off()

pdf("cophenetic.pdf")
coph <- res$measures$cophenetic
plot(2:7,coph, type="b", col="purple")
dev.off()


#Repeat NMF with the selected rank
#Selected rank: 3
res3 <- nmf(tumor_exp,3,nrun=10)





pdf("coefmap.pdf")
coefmap(res3)
dev.off()

pdf("consensusmap.pdf")
consensusmap(res3)
dev.off()

#Compare available NMF methods; brunet is used here
res.multi.method <- nmf(tumor_exp, 3, list('brunet', 'lee', 'ns'), seed=222)
compare(res.multi.method)


W <- basis(res3)
H <- coef(res3)



#Select genes
rank <- 3
seed <- 20220311
tumor_exp=as.data.frame(dat)
#as.data.frame(ing)
mut.nmf <- nmf(tumor_exp, 
               rank = rank, 
               seed = seed,
               nrun=10,
               method = "ns") 

group <- predict(mut.nmf) # Extract predicted subtypes
table(group)
ac=data.frame(group=paste0('NMF_',group))
rownames(ac)=colnames(tumor_exp)
write.table(ac,"group_infor_before_filtering.txt",quote=F,row.names=T)

pdf("1consensusmap_before_filtering.pdf")
consensusmap(mut.nmf)
dev.off()



#Filter genes###
index <- extractFeatures(mut.nmf,"max") 
sig.order <- unlist(index)

length(sig.order)
sig.order<-sig.order[-29]


tumor_exp2 <- tumor_exp[sig.order,]
write.table(rownames(tumor_exp2),"filtered.geneList.txt",quote=F,row.names=F,col.names=F)
write.table(tumor_exp2,"APA_exp_filtered.txt",quote=F)

library(pheatmap)
pheatmap(tumor_exp2,cluster_rows = T,cluster_cols = F)

#Repeat NMF with the selected genes
mut.nmf2 <- nmf(tumor_exp2, 
                rank = rank, 
                seed = seed,
                nrun=10,				
                method = "ns")

group <- predict(mut.nmf2) # Extract predicted subtypes
table(group)

ac=data.frame(group=paste0('NMF_',group))
rownames(ac)=colnames(tumor_exp2)
write.table(ac,"group_infor.txt",quote=F,row.names=T)

pdf("consensusmap.pdf")
consensusmap(mut.nmf2)
dev.off()



#Visualization functions
sample.order <- names(group[order(group)])
#Set colors
jco <- c("#2874C5","#EABF00","#C6524A","#868686")

pdf(file = "consensusmap.pdf")
consensusmap(mut.nmf2,
             labRow = NA,
             labCol = NA,
             annCol = data.frame("cluster"=group[colnames(tumor_exp)]),
             annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
dev.off()

pdf(file = "basismap.pdf" ) 
# Inspect subtype-associated signatures and compare them with the NMF heatmap below

basismap(mut.nmf2,
         cexCol = 1,
         cexRow = 0.3,
         annColors=list(c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])))
dev.off()


group <- predict(mut.nmf) # Extract predicted subtypes
table(group)
sample.order <- names(group[order(group)])

#Normalize the matrix
tumor_exp_normalization<-t(scale(t(tumor_exp)))

tumor_exp_normalization[tumor_exp_normalization>1]=1
tumor_exp_normalization[tumor_exp_normalization<(-1)]=-1


aheatmap(as.matrix(tumor_exp_normalization[,sample.order]), 
         #Rowv=NA, 
         Colv=NA, 
         annCol = data.frame("cluster"=group[sample.order]),
         annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])),
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         revC=TRUE, 
         #cexCol = 0.3, 
         #cexRow = 0.3,
         filename = "NMF_heatmap_before_filtering.pdf")



##################Heatmap of the filtered results####################
group <- predict(mut.nmf2) # Extract predicted subtypes
table(group)
sample.order <- names(group[order(group)])


#Normalize the matrix
tumor_exp_normalization2<-t(scale(t(tumor_exp2)))

tumor_exp_normalization2[tumor_exp_normalization2>1]=1
tumor_exp_normalization2[tumor_exp_normalization2<(-1)]=-1



aheatmap(as.matrix(tumor_exp_normalization2[,sample.order]), 
         #Rowv=NA, 
         Colv=NA, 
         annCol = data.frame("cluster"=group[sample.order]),
         annColors = list(cluster=c("1"=jco[1],"2"=jco[2],"3"=jco[3],"4"=jco[4])),
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         revC=TRUE, 
         #cexCol = 0.3, 
         #cexRow = 0.3,
         filename = "NMF_heatmap.pdf")





ac=data.frame(group=paste0('NMF_',group))
rownames(ac)=colnames(tumor_exp2)
write.table(ac,"group_infor.txt",quote=F,row.names=T)

pheatmap(tumor_exp2,show_colnames = F,annotation_col = ac)
pheatmap(tumor_exp2[,sample.order],
         show_colnames = F,cluster_cols = F,annotation_col = ac)

pheatmap(tumor_exp2,show_colnames = F,annotation_col = ac,
         filename = "NMF_heatmap1.pdf")
pheatmap(tumor_exp2[,sample.order],
         show_colnames = F,cluster_cols = F,annotation_col = ac,
         filename = "NMF_heatmap2.pdf")
dev.off()
