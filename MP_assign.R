library(scalop)
library(tidyverse)
library(Seurat)
library(data.table)
outpath <- "results/"
Mali50_seu.sample.list <- readRDS("Mali50/Mali50_seu.sample.list.rds")
names(Mali50_seu.sample.list)

m_loged_list <- list()
for (i in names(Mali50_seu.sample.list)){ 
  seu <- Mali50_seu.sample.list[[i]]
  seu <- NormalizeData(seu)
  m_loged <- LayerData(seu, assay = "RNA", layer = "data")
  m_loged <- as.matrix(m_loged)
  m_loged_list[[i]] <- m_loged
}

metaprograms_gene <-  readRDS("Mali50/metaprograms_gene.rds")
score_mat <- lapply(names(m_loged_list), function(i){
  m_loged<-m_loged_list[[i]]
  signatures <- scalop::sigScores(m_loged, metaprograms_gene, expr.center = TRUE, conserved.genes = 0.5)
  signatures <-  signatures %>% rownames_to_column(var = "cellid")
  return(signatures)
})
names(score_mat) <- names(m_loged_list)
saveRDS(score_mat,paste0(outpath,"score_mat_list.rds"))

score_df_list <- list()
maxcol_meta_list <- list()
for (i in names(score_mat)){
  score_df <- as.data.frame(score_mat[[i]])
  score_df <-column_to_rownames(score_df, 'cellid')
  maxcol_meta<-maxcol_strict(score_df)
  maxcol_meta<-stack(maxcol_meta)
  setnames(maxcol_meta,2,"MP_assign")
  setnames(maxcol_meta,1,"CellID")
  score_df_list[[i]] <-  score_df
  maxcol_meta_list[[i]] <-  maxcol_meta
}

saveRDS(maxcol_meta_list,paste0(outpath,"maxcol_meta_list.rds"))
maxcol_meta_list[[1]] %>% head()



