library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(dplyr)

combined_seura@meta.data$Classification[combined_seura@meta.data$Classification == ""] <- "Unknown"

unique(combined_seura@meta.data$Classification)
se <- NormalizeData(combined_seura)
se <- FindVariableFeatures(se)
se <- ScaleData(se)
se <- RunPCA(se, npcs = 30)  # adjust number of PCs as needed

## run harmony with default parameters
pbmc <- se %>% RunHarmony("Classification")
pbmc <- pbmc %>% 
  RunHarmony("Classification", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)

p1 <- DimPlot(object = pbmc, reduction = "harmony", pt.size = .1, group.by = "Classification")
p2 <- VlnPlot(object = pbmc, features = "harmony_1", group.by = "Classification",  pt.size = .1)
plot_grid(p1,p2)

DimHeatmap(object = pbmc, reduction = "harmony", cells = 500, dims = 1:3)
pbmc <- pbmc %>%
  FindNeighbors(reduction = "harmony") %>%
  FindClusters(resolution = 0.5) 


#######################

pbmc <- pbmc %>%
  RunTSNE(reduction = "harmony")


p1 <- DimPlot(pbmc, reduction = "tsne", group.by = "Classification", pt.size = .1)
p2 <- DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = .1)
plot_grid(p1, p2)
pbmc <- pbmc %>%
  RunUMAP(reduction = "harmony",  dims = 1:20)

p1 <- DimPlot(pbmc, reduction = "umap", group.by = "Classification", pt.size = .1)
p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE,  pt.size = .1)
plot_grid(p1, p2)

