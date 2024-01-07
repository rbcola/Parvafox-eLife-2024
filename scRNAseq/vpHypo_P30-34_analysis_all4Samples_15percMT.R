### Code derived from https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html 
### Reto B. Cola

rm(list=ls())

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)

set.seed(200591)
outputFolder <- "Z:/data/r.cola/scRNAseq/15percentMT_analysis - Copy/OutputAll4/"

# Load the PBMC dataset
male1.data <- Read10X(data.dir = "Z:/data/r.cola/scRNAseq/GSM4404135_AJ18003_Male1")
male2.data <- Read10X(data.dir = "Z:/data/r.cola/scRNAseq/GSM4404138_AJ19002_Male2")
female1.data <- Read10X(data.dir = "Z:/data/r.cola/scRNAseq/GSM4404136_AJ18004_Female1")
female2.data <- Read10X(data.dir = "Z:/data/r.cola/scRNAseq/GSM4404137_AJ19001_Female2")
# Initialize the Seurat object with the raw (non-normalized data).
male1 <- CreateSeuratObject(counts = male1.data, project = "m1", min.cells = 3, min.features = 200)
male2 <- CreateSeuratObject(counts = male2.data, project = "m2", min.cells = 3, min.features = 200)
female1 <- CreateSeuratObject(counts = female1.data, project = "f1", min.cells = 3, min.features = 200)
female2 <- CreateSeuratObject(counts = female2.data, project = "f2", min.cells = 3, min.features = 200)

#combine the 4 datasets
vpHypo.big <- merge(male1, y = c(male2, female1, female2), add.cell.ids = c("m1", "m2", "f1","f2"), project = "vpHypo")

# Lets examine a few genes in the first thirty cells
vpHypo.big[c("Cck", "Foxb1", "Adcyap1r1", "Adcyap1"), 1:30]
table(vpHypo.big$orig.ident)
GetAssayData(vpHypo.big)[1:10,1:15]

### Standard  pre-processing workflow
## Quality control and cell selection for further analysis

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# PercentageFeatureSet is used to calculate percentage of mitochondiral genes (dying cells exhibit extensive mitochondrial contamination)
vpHypo.big[["percent.mt"]] <- PercentageFeatureSet(vpHypo.big, pattern = "^mt-")

#...and calculate the number of hemoglobin transcripts per cell
Hb_features <- grep(pattern = "^Hb[^(p,e,s)]", x = rownames(x = vpHypo.big[["RNA"]]), value = TRUE)
vpHypo.big[["nHbFeature_RNA"]] <- colSums(x = GetAssayData(object = vpHypo.big, assay = "RNA", slot = "counts")[Hb_features, , drop = FALSE])

# Show QC metrics for the first 5 cells
head(vpHypo.big@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(vpHypo.big, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "nHbFeature_RNA"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(vpHypo.big, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)
plot2 <- FeatureScatter(vpHypo.big, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)
plot1 + plot2 # Patchwork package required for this operation

#based on the plots, we will filter as follows:
vpHypo.big <- subset(vpHypo.big, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 15 & nHbFeature_RNA < 50)
GetAssayData(vpHypo.big)[1:10,1:15]

#Normalize datasets individually by SCTransform(), instead of NormalizeData() prior to integration
vpHypo.list <- SplitObject(vpHypo.big, split.by = "orig.ident")
vpHypo.list <- lapply(X = vpHypo.list, FUN = SCTransform) #conserve.memory helps avoid running out of RAM

features <- SelectIntegrationFeatures(object.list = vpHypo.list, nfeatures = 3000)
vpHypo.list <- PrepSCTIntegration(object.list = vpHypo.list, anchor.features = features)

save.image("milestone_1.RData")

vpHypo.anchors <- FindIntegrationAnchors(object.list = vpHypo.list, normalization.method = "SCT",
                                         anchor.features = features)
save.image("milestone_2.RData")

vpHypo.combined.sct <- IntegrateData(anchorset = vpHypo.anchors, normalization.method = "SCT")

vpHypo.combined.sct <- RunPCA(vpHypo.combined.sct, verbose = FALSE)

ElbowPlot(vpHypo.combined.sct, ndims = 50)

vpHypo.combined.sct <- FindNeighbors(vpHypo.combined.sct, reduction = "pca", dims = 1:30)

vpHypo.combined.sct <- FindClusters(vpHypo.combined.sct, resolution = 0.5)

vpHypo.combined.sct <- RunUMAP(vpHypo.combined.sct, reduction = "pca", dims = 1:30)


save.image("milestone_3.RData")


dimPlot1 <- DimPlot(vpHypo.combined.sct, reduction = "umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)
ggsave(paste0(outputFolder,"dimPlot1.tiff"), plot = dimPlot1, height = 8, width = 8, units = "in")
ggsave(paste0(outputFolder,"dimPlot1.svg"), plot = dimPlot1, height = 8, width = 8, units = "in")
dimPlot2 <- DimPlot(vpHypo.combined.sct, reduction = "umap", group.by = "orig.ident", raster = FALSE)
ggsave(paste0(outputFolder,"dimPlot2.tiff"), plot = dimPlot2, height = 8, width = 8, units = "in")
ggsave(paste0(outputFolder,"dimPlot2.svg"), plot = dimPlot2, height = 8, width = 8, units = "in")

dimHeatMap <- DimHeatmap(vpHypo.combined.sct, dims = 1:5, cells = 500, balanced = TRUE)

# find markers for every cluster compared to all remaining cells, report only the positive ones
vpHypo.combined.sct.markers <- FindAllMarkers(vpHypo.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
vpHypo.combined.sct.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) # note that top_n() has been superseded by slice_min()/slice_max()

saveRDS(vpHypo.combined.sct.markers, "vpHypo_combined_sct_markers.RDS")
#save.image("milestone_4.RData")
#vpHypo.combined.sct.markers <- readRDS("vpHypo_combined_sct_markers.RDS")

# HeatMap for top 20 markers of each cluster
top10 <- vpHypo.combined.sct.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
heatmapPlot <- DoHeatmap(vpHypo.combined.sct, features = top10$gene, raster = FALSE) + NoLegend()
ggsave(paste0(outputFolder,"heatmapPlot.tiff"), plot = heatmapPlot, height = 40, width = 40, units = "in")
#ggsave(paste0(outputFolder,"heatmapPlot.svg"), plot = heatmapPlot, height = 10, width = 10, units = "in")

### Finding differentially expressed features (cluster biomarkers)
# find all markers of your cluster of interest
cluster9.markers <- FindMarkers(vpHypo.combined.sct, ident.1 = 9, min.pct = 0.25)
C9_largestLog2FC <- cluster9.markers %>% arrange(desc(avg_log2FC)) #%>% head(n=20)
saveRDS(C9_largestLog2FC, "C9_largestLog2FC.RDS")

# find all markers distinguishing your cluster of interest from other clusters
cluster9vs5_7_8.markers <- FindMarkers(vpHypo.combined.sct, ident.1 = 9, ident.2 = c(5,7,8), min.pct = 0.25)
cluster9vs5_7_8_largestLog2FC <- cluster9vs5_7_8.markers %>% arrange(desc(avg_log2FC)) #%>% head(n=20)
saveRDS(cluster9vs5_7_8_largestLog2FC, "cluster9vs5_7_8_largestLog2FC.RDS")



featPlot_all <- FeaturePlot(vpHypo.combined.sct, features = c("Cck", "Foxb1","Synpr","Ebf3","Dlk1","Stxbp6"), order = TRUE, raster = FALSE)
ggsave(paste0(outputFolder,"featPlot_all.tiff"), plot = featPlot_all, height = 16, width = 16, units = "in")
ggsave(paste0(outputFolder,"featPlot_all.svg"), plot = featPlot_all, height = 16, width = 16, units = "in")

featPlot_all2 <- FeaturePlot(vpHypo.combined.sct, features = c("Cck", "Foxb1"), blend = TRUE, blend.threshold = 0.25, order=TRUE, raster = FALSE) #+ DarkTheme()
ggsave(paste0(outputFolder,"featPlot_all2.tiff"), plot = featPlot_all2, height = 8, width = 32, units = "in")
ggsave(paste0(outputFolder,"featPlot_all2.svg"), plot = featPlot_all2, height = 8, width = 32, units = "in")

featPlot_PMdCluster <- FeaturePlot(subset(vpHypo.combined.sct, subset= seurat_clusters == 9), features = c("Cck", "Foxb1"), blend = TRUE, blend.threshold = 0.25, order=TRUE, raster = FALSE)
ggsave(paste0(outputFolder,"featPlot_PMdCluster.tiff"), plot = featPlot_PMdCluster, height = 4, width = 16, units = "in")
ggsave(paste0(outputFolder,"featPlot_PMdCluster.svg"), plot = featPlot_PMdCluster, height = 4, width = 16, units = "in")


featPlot_PMd_fixedScale <- FeaturePlot(subset(vpHypo.combined.sct, subset= seurat_clusters == 9), features = c("Cck", "Foxb1"), order = TRUE, raster = FALSE, keep.scale = "all")
ggsave(paste0(outputFolder,"featPlot_PMd_fixedScale.tiff"), plot = featPlot_PMd_fixedScale, height = 4, width = 8, units = "in")
ggsave(paste0(outputFolder,"featPlot_PMd_fixedScale.svg"), plot = featPlot_PMd_fixedScale, height = 4, width = 8, units = "in")

