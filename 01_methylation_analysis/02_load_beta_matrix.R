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

require(GEOquery)
require(Biobase)
library("impute")

# Example group.txt and data.txt contain six methylation arrays in two groups
# These example files illustrate the differential analysis workflow

# Prepare group.txt and data.txt for the GEO dataset of interest
# Alternatively use the GEOquery workflow below

info=read.table("group.txt",sep="\t",header=T)
library(data.table)
b=info
rownames(b)=b[,1]
# If the methylation matrix has already been prepared in a spreadsheet
# Use the fread workflow below
a=fread("data.txt",data.table = F )
a[1:4,1:4]
rownames(a)=a[,1]
a=a[,-1]
beta=as.matrix(a)
beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
a=betaData
a[1:4,1:4]
identical(colnames(a),rownames(b))
# Ensure methylation matrix columns match phenotype rows exactly

library(ChAMP)
# The beta-value matrix must not contain missing values
myLoad=champ.filter(beta = a,pd = b)
myLoad
save(myLoad,file = 'step1-output.Rdata')

# If downloading the methylation matrix from GEO
# Review the following dataset-specific code
if(F){
  require(GEOquery)
  require(Biobase)
  eset <- getGEO("GSE68777",destdir = './',AnnotGPL = T,getGPL = F)
  beta.m <- exprs(eset[[1]])
  ## Prepare clinical metadata; adapt the selected columns for each study
  pD.all <- pData(eset[[1]])
  pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2")]
  head(pD)
  names(pD)[c(3,4)] <- c("group", "sex")
  pD$group <- sub("^diagnosis: ", "", pD$group)
  pD$sex <- sub("^Sex: ", "", pD$sex)
  
  library(ChAMP)
  # The beta-value matrix must not contain missing values
  myLoad=champ.filter(beta = beta.m ,pd = pD)
  myLoad
  save(myLoad,file = 'step1-output.Rdata')
}

# Both approaches create a ChAMP input object
# Use myLoad for subsequent analyses