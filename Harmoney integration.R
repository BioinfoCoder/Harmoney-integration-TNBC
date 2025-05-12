library(Seurat)
library(harmony)
library(patchwork)
# Sample IDs
samples <- c("CID44971", "CID4465", "1160920F", "1142243F")#CID4535 & CID4290 ER

# Base path to parent directory (edit this path!)
raw_dir <- "C:/Users/doaad/Downloads/4739739/raw_count_matrices"
Filt_dir<-"C:/Users/doaad/Downloads/4739739/filtered_count_matrices"
meta_dir<-"C:/Users/doaad/Downloads/4739739/metadata"

seurat_list <- list()

for (sample in samples) {
  message("Processing ", sample, "...")
  
  # Construct file paths
  filt_path <- file.path(Filt_dir, paste0(sample, "_filtered_count_matrix"))
  raw_path  <- file.path(raw_dir, paste0(sample, "_raw_feature_bc_matrix"))
  meta_path <- file.path(meta_dir, paste0(sample, "_metadata.csv"))
  
  # Load data
  filt_counts <- Read10X(data.dir = filt_path, gene.column = 1)
  raw_counts  <- Read10X(data.dir = raw_path)
  metadata    <- read.csv(meta_path, row.names = 1)
  
  # Match metadata to barcodes
  metadata <- metadata[colnames(filt_counts), , drop = FALSE]
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = filt_counts, meta.data = metadata, project = sample)
  
  # Subset and add raw counts as a new assay
  raw_counts <- raw_counts[, colnames(filt_counts)]
  seu[["Raw"]] <- CreateAssayObject(counts = raw_counts)
  
  # Add sample ID
  seu$sample_id <- sample
  
  # ----------- Preprocessing ----------- #
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
  seu <- ScaleData(seu, features = VariableFeatures(seu))
  seu <- RunPCA(seu, features = VariableFeatures(seu))
  
  # UMAP before Harmony
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:20)
  umap_pre <- DimPlot(seu, reduction = "umap", group.by = "sample_id") + ggtitle(paste0(sample, ": UMAP (Before Harmony)"))
  
  # Harmony
  #seu <- RunHarmony(seu, group.by.vars = "sample_id")
  
  # UMAP after Harmony
  # seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)
  #umap_post <- DimPlot(seu, reduction = "umap", group.by = "sample_id") + ggtitle(paste0(sample, ": UMAP (After Harmony)"))
  
  # Combine and print
  plot_combined <- umap_pre
  print(plot_combined)
  
  # Store in list
  seurat_list[[sample]] <- seu
}

# Combine all Seurat objects into one
combined_seura <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)], add.cell.ids = names(seurat_list))

# ----------- Preprocessing ----------- #
# Normalize the data
combined_seura <- NormalizeData(combined_seura)
combined_seura <- FindVariableFeatures(combined_seura, selection.method = "vst", nfeatures = 2000)
combined_seura <- ScaleData(combined_seura, features = VariableFeatures(combined_seura))

# Perform PCA on combined data
combined_seura <- RunPCA(combined_seura, features = VariableFeatures(combined_seura))

# Apply Harmony across all samples using PCA reduction
combined_seura <- RunHarmony(combined_seura, group.by.vars = "sample_id")

# Visualize PCA before and after Harmony in the same plot
# PCA before Harmony
p1 <- DimPlot(combined_seura, reduction = "pca", group.by = "sample_id", label = TRUE) + ggtitle("PCA Before Harmony")

# PCA after Harmony
p2 <- DimPlot(combined_seura, reduction = "harmony", group.by = "sample_id", label = TRUE) + ggtitle("PCA After Harmony")

# Combine the two plots
library(patchwork)
p1 + p2
samples <- c("CID44971", "CID4465", "1160920F", "1142243F")

# Filter, get unique combinations, sort by patient ID, and print
df <- combined_seura@meta.data[
  combined_seura@meta.data$patientid %in% samples, 
  c("patientid", "Classification")
]

df <- unique(df)
df <- df[order(df$patientid), ]
print(df)
