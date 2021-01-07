# Copyright (c) [2021] [Ricardo O. Ramirez Flores, Christian H. Holland]
# roramirezf@uni-heidelberg.de; christian.holland@bioquant.uni-heidelberg.de

library(here)
library(tidyverse)
library(Seurat)

# Manual transformation from Seurat v2 to v3 ------------------------------
seurat_v2 <- get(load(here("data/meta-mouse-vs-human/tissue.rdata")))
rm(tissue)

# extract relevant slots from seurat object version 2
counts <- as.matrix(seurat_v2@raw.data)
meta <- seurat_v2@meta.data
data <- as.matrix(seurat_v2@data)

rm(seurat_v2)

# construct seurat object version 3
seurat_v3 <- CreateSeuratObject(counts = counts, min.cells = 3,
                                min.features = 200)
seurat_v3@assays$RNA@data <- data
seurat_v3@meta.data <- meta


# UMAP analysis -----------------------------------------------------------
DefaultAssay(seurat_v3) = "RNA"

seurat_v3 = seurat_v3 %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA(features = VariableFeatures(.)) %>%
  RunUMAP(dims = 1:10)



# DEA between cirrhotic and healthy patients ------------------------------

get_condition_markers <- function(ct, data, expr, sample_ann) {

  ct_cells <- which(sample_ann$annotation_lineage == ct)

  ct_ann <- sample_ann[ct_cells,]

  red_data <- data[, ct_cells]
  red_expr <- expr[, ct_cells]

  red_obj <- CreateSeuratObject(counts = red_expr,
                                min.cells = 3,
                                min.features = 200)

  red_obj@assays$RNA@data <- red_data
  red_obj@meta.data <- ct_ann

  DefaultAssay(red_obj) <- "RNA"
  Idents(red_obj) <- "condition"

  FindAllMarkers(object = red_obj,
                 assay = "RNA",
                 logfc.threshold = 0.5,
                 only.pos = F)
}

meta <- meta %>%
  rownames_to_column("cell_id")

dea_res <- map(set_names(unique(meta$annotation_lineage)),
               get_condition_markers,
               data = data,
               expr = expr,
               sample_ann = meta)

saveRDS(dea_res, here("data/meta-mouse-vs-human/single_cell_degs.rds"))
