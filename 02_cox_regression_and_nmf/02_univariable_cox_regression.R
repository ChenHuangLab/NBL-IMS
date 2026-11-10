myLoad1=champ.filter(beta =myLoad$beta ,pd = pD) 
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5) 
group_list=myLoad$pd$group 
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)


Coxoutput <- NULL 
for(i in 3:ncol(com)){
  g <- colnames(com)[i]
  cox <- coxph(Surv(OS_MONTHS,OS_STATUS) ~ com[,i], data = com) # Univariable Cox proportional hazards model
  coxSummary = summary(cox)
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1]))
}
pcutoff <- 0.05
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),]

betaFrame#Exclude normal samples
betaF

rame=read.csv('C:/Users/pc/Desktop/NBL_OS_pD.csv', header = TRUE, sep = ",", quote = "\"",
          dec = ".", fill = TRUE, comment.char = "")



Coxoutput <- NULL 
for(i in 3:ncol(h)){
  g <- colnames(h)[i]
  cox <- coxph(Surv(OS_DAYS,OS_STATUS) ~ h[,i], data =h ) # Univariable Cox proportional hazards model
  coxSummary = summary(cox)
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1]))
}
pcutoff <- 0.05
topgene1 <- Coxoutput[which(Coxoutput$pvalue < pcutoff),]


save(cox, Coxoutput, coxSummary, topgene, topgene1,
     file = "10.18cox_main.RData")
