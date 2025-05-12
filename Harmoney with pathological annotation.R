library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 30)  # adjust number of PCs as needed

## run harmony with default parameters
pbmc <- seu %>% RunHarmony("Classification")
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
sessionInfo()
FeaturePlot(object = pbmc, features= c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"), 
            min.cutoff = "q9", cols = c("lightgrey", "blue"), pt.size = 0.5)
