load("D:/RProject/NBL/metha/10.23_immunepair_cpG&gene.Rdata")
load("D:/RProject/NBL/metha/10.18cox_main.RData")
load("D:/RProject/NBL/metha/10.23_immune_paired_sep_pvalue.Rdata")

intersect(topgene1$gene,immunepair$cpG)
