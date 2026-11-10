Coxoutput <- NULL 
for(i in 4:ncol(d)){
  g <- colnames(d)[i]
  cox <- coxph(Surv(OS_DAYS,OS_STATUS) ~ d[,i], data =d ) # Univariable Cox proportional hazards model
  coxSummary = summary(cox)
  Coxoutput <- rbind.data.frame(Coxoutput,
                                data.frame(gene = g,
                                           HR = as.numeric(coxSummary$coefficients[,"exp(coef)"])[1],
                                           z = as.numeric(coxSummary$coefficients[,"z"])[1],
                                           pvalue = as.numeric(coxSummary$coefficients[,"Pr(>|z|)"])[1]))
}
pcutoff <- 0.05
topgene <- Coxoutput[which(Coxoutput$pvalue < pcutoff),]

for(i in 3:ncol(n)){
  g <- colnames(n)[i]
  cox2 <- coxph(Surv(OS_DAYS,OS_STATUS) ~ n[,i], data =n )}
test=test[,-212]
View(test)
test=test[,-202]
test=test[,-183]
test=test[,-143]
test=test[,-108]
test=test[,-87]
test=test[,-79]
test=test[,-75]
test=test[,-43]
test=test[,-40]
test=test[,-29]
test=test[,-14]
test=test[,-12]