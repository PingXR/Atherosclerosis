library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
setwd("~/20240103_Atherosis/result/1-dealdata")

#GEO数据库
data <- read.csv(file="~/Atherosis_0723/data_human/GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt.gz", 
                 sep="\t",header = T,row.names = 1)
seurat.obj_old <- CreateSeuratObject(counts = data,project = "human_AS")

#SRA数据
counts1 <-
  Read10X(data.dir = "/home/pingxr/Atherosis_0723/data_human/FASTQ/cellranger-result/rca1-1/outs/filtered_feature_bc_matrix/")
seurat_obj1 <- CreateSeuratObject(counts = counts1,project = "01") 
seurat_obj1[["group"]] <- "Atherosis"
seurat_obj1[["donor"]] <- "rac1-1"

counts2 <-
  Read10X(data.dir = "/home/pingxr/Atherosis_0723/data_human/FASTQ/cellranger-result/rca1-2/outs/filtered_feature_bc_matrix/")
seurat_obj2 <- CreateSeuratObject(counts = counts2,project = "02") 
seurat_obj2[["group"]] <- "Atherosis"
seurat_obj2[["donor"]] <- "rac1-2"

counts3 <-
  Read10X(data.dir = "/home/pingxr/Atherosis_0723/data_human/FASTQ/cellranger-result/rca2-1/outs/filtered_feature_bc_matrix/")
seurat_obj3 <- CreateSeuratObject(counts = counts3,project = "03") 
seurat_obj3[["group"]] <- "Atherosis"
seurat_obj3[["donor"]] <- "rac2-1"

counts4 <-
  Read10X(data.dir = "/home/pingxr/Atherosis_0723/data_human/FASTQ/cellranger-result/rca2-2/outs/filtered_feature_bc_matrix/")
seurat_obj4 <- CreateSeuratObject(counts = counts4,project = "04") 
seurat_obj4[["group"]] <- "Atherosis"
seurat_obj4[["donor"]] <- "rac2-2"

counts5 <-
  Read10X(data.dir = "/home/pingxr/Atherosis_0723/data_human/FASTQ/cellranger-result/rca3-1/outs/filtered_feature_bc_matrix/")
seurat_obj5 <- CreateSeuratObject(counts = counts5,project = "05") 
seurat_obj5[["group"]] <- "Atherosis"
seurat_obj5[["donor"]] <- "rac3-1"

counts6 <-
  Read10X(data.dir = "/home/pingxr/Atherosis_0723/data_human/FASTQ/cellranger-result/rca3-2/outs/filtered_feature_bc_matrix/")
seurat_obj6 <- CreateSeuratObject(counts = counts6,project = "06") 
seurat_obj6[["group"]] <- "Atherosis"
seurat_obj6[["donor"]] <- "rac3-2"

counts7 <-
  Read10X(data.dir = "/home/pingxr/Atherosis_0723/data_human/FASTQ/cellranger-result/rca3-3/outs/filtered_feature_bc_matrix/")
seurat_obj7 <- CreateSeuratObject(counts = counts7,project = "07") 
seurat_obj7[["group"]] <- "Atherosis"
seurat_obj7[["donor"]] <- "rac3-3"

counts8 <-
  Read10X(data.dir = "/home/pingxr/Atherosis_0723/data_human/FASTQ/cellranger-result/rca4-1/outs/filtered_feature_bc_matrix/")
seurat_obj8<- CreateSeuratObject(counts = counts8,project = "08") 
seurat_obj8[["group"]] <- "Atherosis"
seurat_obj8[["donor"]] <- "rac4-1"
# 合并8个样本
merged <- merge(seurat_obj1,c(seurat_obj2,seurat_obj3,seurat_obj4,
                              seurat_obj5,seurat_obj6,seurat_obj7,seurat_obj8))
saveRDS(merged,"merged_SRA.rds")

# 取交
#GEO数据库
data <- read.csv(file="~/Atherosis_0723/data_human/GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt.gz", 
                 sep="\t",header = T,row.names = 1)
seurat.obj_GEO <- CreateSeuratObject(counts = data,project = "human_AS")
#SRA原始数据
seurat.obj_SRA<- readRDS("merged_SRA.rds")
seurat.obj_SRA[["barcode"]] <- colnames(seurat.obj_SRA)
seurat.obj_GEO[["barcode"]] <- colnames(seurat.obj_GEO)
library(stringr)
seurat.obj_SRA$barcode <- str_replace_all(seurat.obj_SRA$barcode,"-1_",".")
seurat.obj <- subset(seurat.obj_SRA,subset = barcode %in% seurat.obj_GEO$barcode)
ct <- seurat.obj@assays$RNA@counts
sce2 <-  CreateSeuratObject(counts = ct,min.cells = 6)
sce2@meta.data$orig.ident <- seurat.obj@meta.data$orig.ident
sce2[['group']] <- seurat.obj@meta.data$group
sce2[['donor']] <- seurat.obj@meta.data$donor
saveRDS(sce2,"rawseurat.rds")



