suppressPackageStartupMessages({
  library(arrow)
  library(Seurat)
  library(spacexr)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

st_data_path <- args[1]
ref_data_file <- args[2]
output_path <- args[3]

# Build reference object-------------------------------------------------------------------------------------------------------
sc_ref_data <- readRDS(ref_data_file)
counts <- Seurat::GetAssayData(sc_ref_data,assay = "RNA",slot = "counts")
cluster <- as.factor(sc_ref_data$sub_cell_type_final)
cluster <- droplevels(cluster)
names(cluster) <- colnames(sc_ref_data)
nUMI <- sc_ref_data$nCount_RNA
names(nUMI) <- colnames(sc_ref_data)
set.seed(123) 
reference <- Reference(counts, cluster, nUMI)

# query --------------------------------------------------------------------------------------------
st_data <- readRDS(st_data_path)
dim(st_data)

counts_st <- st_data[["Spatial"]]$counts
cells_st <- colnames(st_data[["Spatial"]])
coords <- Seurat::GetTissueCoordinates(object = st_data)
coords <- coords[cells_st,c("x","y")]
nUMI <- colSums(counts_st)
query <- SpatialRNA(coords, counts_st, nUMI)

# RCTD --------------------------------------------------------------------------------------------
myRCTD <- create.RCTD(query, reference, max_cores = 64)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
saveRDS(myRCTD, file = paste0(output_path, "rctd.rds"))
