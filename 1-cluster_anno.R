setwd("~/20240103_Atherosis/v2/result/1-dealdata/figure")
#### 导入包 ----
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
library(scater)
library(hdf5r)

path <- "~/20240103_Atherosis/result/1-dealdata"
seurat.obj <- readRDS(paste0(path,"/rawseurat.rds"))
# 质控----
seurat.obj <- PercentageFeatureSet(seurat.obj,
                                   pattern = "^MT-",
                                   col.name = "percent.mt")
Idents(seurat.obj) <- "donor"
VlnPlot(
  seurat.obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  stack = TRUE,
  pt.size = 0,
  flip = TRUE,
  col = c("#A5D38F", "#E45D61", "#96C3D8" )
)
seurat.obj <- subset(seurat.obj, subset =
                       nFeature_RNA >= 200 &
                       nFeature_RNA <= 3500 &  
                       percent.mt <= 7)    
dim(seurat.obj)
VlnPlot(
  seurat.obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  stack = TRUE,
  pt.size = 0,
  flip = TRUE,
  col = c("#A5D38F", "#E45D61", "#96C3D8" )
)
seurat.obj <- NormalizeData(seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
seurat.obj <- ScaleData(seurat.obj)
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))

#### 降维聚类分群 ----
ElbowPlot(seurat.obj, ndims = 50)
seurat.obj <- RunTSNE(seurat.obj, dims = 1:50)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:50)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:50)
seurat.obj <- FindClusters(seurat.obj,
                           resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
library(clustree)
clustree(seurat.obj@meta.data,prefix = "RNA_snn_res.")
ggsave("clustree.pdf",width=8,height=10)
seurat.obj <- FindClusters(seurat.obj,
                           resolution = 0.5)
DimPlot(seurat.obj,
        reduction = "umap",
        label = T,
        pt.size = 0.5)
ggsave("umap.pdf",width=6,height=5)
DimPlot(seurat.obj,
        reduction = "tsne", # tsne, umap, pca
        label = T,
        pt.size = 0.5)
ggsave("tsne.pdf",width=6,height=5)


#### marker ----
DefaultAssay(seurat.obj) <- "RNA"
seurat.obj <- NormalizeData(seurat.obj)
features1 <- c(
  "ACTA2","TAGLN","CNN1","TNS1",
  "C7","LUM","DCN",
  "SERPINE2","OMD",
  "CD163","CD14","RNASE1",
  "CCL21","CD36",
  "RERGL","NET1",
  "PECAM1","CLDN5",
  "MS4A1","CD79B", 
  "GAS6","TNFRSF11B",
  "CD3D","CD3E","TRAC",
  "TXNDC5","EAF2",
  "IGHG1","IGHG2",
  "KLRD1","PRF1","GNLY", 
  "S100B","PLP1","GPM6B",
  "TPSAB1","CPA3"
  # "COL1A1", "PDGFRA", "BSG","COL3A1","GPC3", 
  # "PECAM1", "CD31", "CLDN5", 
  # "CD14", "CD163", "MARCO", "MSR1", "MRC1", 
  # "CD3D","CD3E","CD3G", 
  # "CNN1", "ACTA2", "TAGLN", "RGS5", "DES", "LGR6",  
  # "CSPG4", "TRPC6", "PDGFRB", 
  # "CD79A", "CD24", "MS4A1", "CD19", 
  # "CD27", "SLAMF7", 
  # "MS4A2", "CPA3", "KIT","TPSAB1",  
  # "KLRD1","PRF1","GNLY","NKG7","NCAM1","FCGR3A", 
  # "TUBB3","SYN1" 
)
FeaturePlot(
  seurat.obj,
  features = features1,
  reduction = "tsne",
  order = T,
  ncol = 4,
  min.cutoff = 0
)
ggsave("featureplot2.pdf",width=24,height=30)
DotPlot(seurat.obj, features = features1, col.min = 0)+ theme(axis.text.x = element_text(angle = 45,vjust = 0))
ggsave("dotplot2.pdf",width=10,height=6)
markers <- FindAllMarkers(
  seurat.obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers,"findallmarkers_0.5.csv")

saveRDS(seurat.obj,"~/20240103_Atherosis/v2/result/1-dealdata/seurat_unanno.rds")

#### 细胞注释 ----
seurat.obj <- readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_unanno.rds")
# 分配细胞名字
new.cluster.ids <- c(
  "Fibroblast 1", # 00  Fibroblast
  "Endothelial",# 01   Endothelial
  "Macrophage", # 02  Macrophage
  "Fibromyocyte", # 03  Fibromyocyte
  "T cell", # 04  T cell
  "Smooth muscle cell",# 05 
  "Pericyte 1",# 06 Pericyte 1
  "Pericyte 2", # 07 Pericyte 2 
  "B cell",#08 B cell
  "C9_APOE",#09
  "C10_S100A8",#10
  "Plasma cell 1",#11 Plasma cell 1
  "C12_EDN1",#12
  "Fibroblast 2",#13 Fibroblast 2
  "Neuron",#14
  "Plasma cell 2",#15 Plasma cell 2
  "NK cell",#16 NK cell
  "C17_GZMB",#17
  "Mast cell",#18 Mast cell
  "C19_TFF3"#19
)
names(new.cluster.ids) <- levels(seurat.obj)
seurat.obj <- RenameIdents(seurat.obj,
                           new.cluster.ids)
seurat.obj$cell_type <- Idents(seurat.obj)
DimPlot(seurat.obj,
        label = T,
        repel = T,
        pt.size = 0.3,
        reduction = "tsne") 
DimPlot(seurat.obj,
        label = T,
        repel = T,
        pt.size = 0.5,
        reduction = "umap") 

#### 保存结果
setwd("~/20240103_Atherosis/v2/result/1-dealdata")
saveRDS(seurat.obj,file = "seurat_integration_anno2.rds")
head(seurat.obj)
table(seurat.obj$seurat_clusters)
table(seurat.obj$cell_type)
table(seurat.obj$donor)
seurat.obj
dim(seurat.obj)







