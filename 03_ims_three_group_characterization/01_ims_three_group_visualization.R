##Input notes
##cluster contains group assignments
##geneID contains unique identifiers; retain it for differential analysis


rm(list=ls())
??pacman

# Install pacman
#install.packages("pacman")
# Load pacman
library(pacman)
# Load packages with pacman
#install.packages("AnnotationDbi")
#library(AnnotationDbi)
pacman::p_load(tidyverse, ggrepel, reshape2, rio)
# Load data
#DEG_limma_voom <- read.table("DEG_limma_voom.txt", header=T, sep="\t")
load("D:/RProject/NBL/catecholamine_related/three_group_markers/NMF_3grp.matrix.rds")

DEG_limma_voom=NMF_3grp

# Filter nonsignificant results
d1 <- 
  DEG_limma_voom %>%  
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No change'))) %>% 
  filter(change != 'No change') %>% 
  select(-change)

# Create volcano-plot groups
set.seed(123)
g <- 
  rep(paste0('Volcano_', 1:10),length.out = nrow(d1)) %>% .[sample(nrow(d1))]

g=NMF_3grp[rownames(d1),]
g=g$cluster


g[1:50]

# Define P-value groups
data <- 
  d1 %>% 
  rownames_to_column('gene') %>% 
  mutate(group = g) %>% 
  mutate(label = ifelse(adj.P.Val < 0.05,"adjust P-val<0.05","adjust P-val >= 0.05"))

table(data$group)
table(data$label)




#Select the top 10 genes per group by absolute log fold change
TopGene <- 
  data %>% 
  group_by(group) %>% 
  distinct(gene, .keep_all = T) %>% 
  top_n(10, abs(logFC))
table(TopGene$group)
TopGene$geneID



# Prepare background bar data
dbar <- 
  data %>% 
  group_by(group) %>% 
  summarise_all(list(min = min, max = max)) %>% 
  select(group, logFC_min, logFC_max, label_min) %>% 
  rename(label = label_min)
dbar



# Generate the plot
ggplot()+
  geom_col(data = dbar,  # Draw negative background bars
           mapping = aes(x = group,y = logFC_min),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_col(data = dbar, # Draw positive background bars
           mapping = aes(x = group,y = logFC_max),
           fill = "#dcdcdc",alpha = 0.6, width = 0.7) +
  geom_jitter(data = data, # Draw all data points
              aes(x = group, y = logFC, color = label),
              size = 0.85,
              width =0.3) +
  geom_jitter(data = TopGene, # Highlight the top 10 data points
              aes(x = group, y = logFC, color = label),
              size = 1,
              width =0.35) +
  geom_tile(data = TopGene, # Draw central group annotations
            aes(x = group,
                y = 0,
                fill = group),
            height=1.5,
            color = "black",
            alpha = 0.6,
            show.legend = F) +
  ggsci::scale_fill_npg() + # Set custom colors
  ggsci::scale_color_npg() + # Set custom colors
  geom_text_repel(data = filter(data, geneID %in% TopGene$geneID),  # Filter genes to label
                  aes(x = group, y = logFC, label = geneID),
                  size = 2, 
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 15),
                  color = 'black',
                  force = 1.2,
                  arrow = arrow(length = unit(0.008, "npc"),
                                type = "open", ends = "last")) +
  labs(x="Cluster", y="Average logFC") +
  geom_text(data=TopGene, # Draw central group labels
            aes(x=group, 
                y=0, 
                label=group),
            size = 3,
            color ="white") +
  theme_minimal() + 
  theme(axis.title = element_text(size = 13,color = "black",face = "bold"),
        axis.line.y = element_line(color = "black",size = 1.2),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 13))
ggsave('pic.png', width = 4, height = 6, bg = 'white')
ggsave('3nmf_deg.pdf', width = 4, height = 6, bg = 'white')
