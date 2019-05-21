# https://github.com/satijalab/seurat/issues/1020

set.seed(1234567)
library(Seurat)
library(dplyr)
library(cowplot)

pbmc.data <- Read10X(data.dir = "/mnt/scratchb/bioinformatics/baller01/20190121_BallereauS_BC_ScRsCrs/filtered_gene_bc_matrices/hg19/")

pbmc <- CreateSeuratObject(counts = pbmc.data, min.cells = 3, min.features  = 200, project = "10X_PBMC", assay = "RNA")

mito.genes <- grep(pattern = "^MT-", x = rownames(pbmc@assays[["RNA"]]), value = TRUE)

percent.mito <- Matrix::colSums(pbmc@assays[["RNA"]][mito.genes, ])/Matrix::colSums(pbmc@assays[["RNA"]])

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats

# Seurat v2 function, but shows compatibility in Seurat v3
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito") 
# in case the above function does not work simply do:
pbmc$percent.mito <- percent.mito

# [v3]
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# We filter out cells that have unique gene counts (nFeature_RNA) over 2,500 or less than
# 200 Note that > and < are used to define a 'gate'.  
# -Inf and Inf should be used if you don't want a lower or upper threshold.
pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mito >  -Inf & percent.mito < 0.05 )
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
head(x = HVFInfo(object = pbmc))

 # or with VariableFeatures()
 # Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

 # plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nCounts_RNA", "percent.mito"))
pbmc <- RunPCA(object = pbmc,  npcs = 30, verbose = FALSE)
# Examine and visualize PCA results a few different ways
DimPlot(object = pbmc, reduction = "pca")
# Dimensional reduction plot, with cells colored by a quantitative feature
FeaturePlot(object = pbmc, features = "MS4A1")
# Scatter plot across single cells, replaces GenePlot
FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "PC_1")
FeatureScatter(object = pbmc, feature1 = "MS4A1", feature2 = "CD3D")
# Scatter plot across individual features, repleaces CellPlot
CellScatter(object = pbmc, cell1 = "AGTCTACTAGGGTG", cell2 = "CACAGATGGTTTCT")
VariableFeaturePlot(object = pbmc)
# Violin and Ridge plots
VlnPlot(object = pbmc, features = c("LYZ", "CCL5", "IL32"))
RidgePlot(object = pbmc, feature = c("LYZ", "CCL5", "IL32"))
# Heatmaps
DimHeatmap(object = pbmc, reduction = "pca", cells = 200, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for
# expediency.  More approximate techniques such as those implemented in
# PCElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(object = pbmc, reduction = "pca", dims = 20, num.replicate = 100,  prop.freq = 0.1, verbose = FALSE)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20, reduction = "pca")
JackStrawPlot(object = pbmc, dims = 1:20, reduction = "pca")
ElbowPlot(object = pbmc)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5, algorithm = 1)
head(pbmc@meta.data)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE)
# note that you can set do.label=T to help label individual clusters
DimPlot(object = pbmc, reduction = "tsne")
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
DimPlot(pbmc, reduction = "umap", split.by = "seurat_clusters")

