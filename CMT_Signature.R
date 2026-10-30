library(tidyverse)
library(data.table)

CM_DEGs <- read.csv("DESeq2_all_res.csv")
head(CM_DEGs)

top_DEGs_filtered <- CM_DEGs %>%
  filter(padj < 0.05) %>%      
  filter(log2FoldChange > 1.25) %>%     
  group_by(CMT) %>%                  
  arrange(padj) %>%                  
  slice_head(n = 75) %>%             
  ungroup()
top_DEGs_filtered

saveRDS(top_DEGs_filtered, "top_DEGs_filtered.rds")