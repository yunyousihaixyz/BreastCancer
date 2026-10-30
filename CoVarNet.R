library(CoVarNet)
library(data.table)
library(tidyverse)
meta <- readRDS("meta.rds")

rt_sp <- names(table(meta$sampleID))[table(meta$sampleID) >= 100]
meta <- meta[meta$sampleID %in% rt_sp, ]
head(meta)

mat_fq_raw<-freq_calculate(meta)
mat_fq_norm<-freq_normalize(mat_fq_raw,normalize="minmax")
dim(mat_fq_norm)

res <- nmf(mat_fq_norm, rank = 2:20, method = "nsNMF", seed = rep(123456, 6), .options = "vp")
plot(res)

K=5
NMF_K<-nmf(mat_fq_norm, K, method = "nsNMF", seed = rep(77, 6), nrun = 30)
NMF_K

colnames(basis(NMF_K))=c("CM01","CM02","CM03","CM04","CM05")
rownames(coef(NMF_K))=c("CM01","CM02","CM03","CM04","CM05")

gr.weight_all(NMF_K)
gr.weight_top(NMF_K,num=15)
