# recluster----
library(Seurat)
setwd("~/20240103_Atherosis/result/Fig_SMC/subset_seurat")

seurat.obj <- readRDS("~/20240103_Atherosis/result/1-dealdata/seurat_integration_anno.rds")
sub.seurat.obj <- subset(
  seurat.obj,
  idents = "SMC"
)
saveRDS(sub.seurat.obj, file = "SMC亚群.rds")
DefaultAssay(sub.seurat.obj) <- "RNA"
sub.seurat.obj <- FindVariableFeatures(sub.seurat.obj,
                                       selection.method = "vst",
                                       nfeatures = 2000)
all.genes <- rownames(sub.seurat.obj)
sub.seurat.obj <- ScaleData(sub.seurat.obj, features = all.genes)
sub.seurat.obj <- RunPCA(sub.seurat.obj,
                         features = VariableFeatures(sub.seurat.obj))
ElbowPlot(sub.seurat.obj, ndims = 50)
ndims <- 15
sub.seurat.obj <- FindNeighbors(sub.seurat.obj, dims = 1:ndims)
sub.seurat.obj <- FindClusters(sub.seurat.obj, resolution = 0.3)
sub.seurat.obj <- RunUMAP(sub.seurat.obj, dims = 1:ndims)
sub.seurat.obj <- RunTSNE(sub.seurat.obj, dims = 1:ndims)
DimPlot(sub.seurat.obj, reduction = "umap")
DefaultAssay(sub.seurat.obj) <- "RNA"
FeaturePlot(sub.seurat.obj,
            features = c(
              "ACTA2", "PHACTR1", "PDGFR",  
              "TNFRSF11B",
              "SOX9",
              "C3",
              "CD68", "LGALS3"
            ))
ggsave("marker.pdf",height = 10,width = 10)
markers <- FindAllMarkers(
  sub.seurat.obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
top.markers <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)
write.table(top.markers[top.markers$cluster=="0",]$gene,
            quote = F, row.names = F)

DotPlot(sub.seurat.obj, features = unique(top.markers$gene),
        col.min = 0) +
  theme(axis.text.x = element_text(size = 12,angle = 45,hjust = 1) )
ggsave("dotplot.pdf",height = 5,
       width = 10) 
new.cluster.ids <- c(
  "SMC1", 
  "SMC2"
)
names(new.cluster.ids) <- levels(sub.seurat.obj)
sub.seurat.obj <- RenameIdents(sub.seurat.obj, new.cluster.ids)
sub.seurat.obj$sub_cell_type <- Idents(sub.seurat.obj)

DimPlot(sub.seurat.obj, reduction = "umap", label = T) + NoLegend()
head(sub.seurat.obj)
saveRDS(sub.seurat.obj, 
        file = "SMC_anno.rds")

# tsne----
library(Seurat)
library(tidyverse)
setwd("~/20240103_Atherosis/result/Fig_SMC/subset_seurat")
path <-
  "SMC_anno.rds"
seurat_obj <- readRDS(path)
DefaultAssay(seurat_obj) <- "RNA"
tsne = seurat_obj@reductions$tsne@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = seurat_obj$sub_cell_type)
head(tsne)
a <- DimPlot(
  seurat_obj,
  reduction = "tsne",
  group.by = "sub_cell_type",
  label = F,
  pt.size = 2.5,
  label.size = 4,
  raster = TRUE
) + 
  scale_color_manual(values = c("#e64b35ff","#4dbbd5ff","#6B70B0","#3D6AAA","#394D9B","#75ACC3","#20ACBD","#38509F"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"),
  )+
  theme(
    legend.title = element_blank(), 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=10), 
  ) +  
  guides(color = guide_legend(override.aes = list(size=4)))+ 
  geom_segment(aes(x = min(tsne$tSNE_1) , y = min(tsne$tSNE_2) ,
                   xend = min(tsne$tSNE_1) +10, yend = min(tsne$tSNE_2) ),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.1,"inches")))+ 
  geom_segment(aes(x = min(tsne$tSNE_1)  , y = min(tsne$tSNE_2)  ,
                   xend = min(tsne$tSNE_1) , yend = min(tsne$tSNE_2) + 10),
               colour = "black", size=0.5,arrow = arrow(length = unit(0.1,"inches"))) +
  annotate("text", x = min(tsne$tSNE_1) +3, y = min(tsne$tSNE_2) -1.5, label = "tsne_1",
           color="black",size = 2, ) + 
  annotate("text", x = min(tsne$tSNE_1) -1.5, y = min(tsne$tSNE_2) + 3, label = "tsne_2",
           color="black",size = 2, angle=90) 
ggsave(
  file.path( "cell_type_tsne2.pdf"),
  plot = a,
  height =4,
  width = 5
)

# monocle2----
library(monocle)
library(Seurat)
library(dplyr) 
library(ggsci)
setwd("~/20240103_Atherosis/result/Fig_SMC/Supply")
seurat.obj <-
  readRDS("~/20240103_Atherosis/result/Fig_SMC/subset_seurat/SMC_anno.rds")
table(seurat.obj@meta.data$sub_cell_type)
DefaultAssay(seurat.obj) <- 'RNA'
data <- as(as.matrix(seurat.obj@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat.obj@meta.data)
fData <-
  data.frame(gene_short_name = row.names(data),
             row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
data <- data[,rownames(pd)]
cds <- newCellDataSet(
  data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <-
  detectGenes(cds, min_expr = 0.1) 
head(fData(cds))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) 
deg <- FindAllMarkers(seurat.obj,
                      assay = "RNA",only.pos = T)
top_deg <- deg %>% 
  filter(p_val_adj < 0.25) %>%  
  arrange("avg_log2FC") %>% 
  group_by(cluster) %>% 
  top_n(50, wt = avg_log2FC) %>%   
  pull(gene)
write.table(
  deg,
  file = "SMC_monocle_deg.txt",
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F
)
ordergene <- top_deg
cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = 2)
saveRDS(cds, "SMC_cds.rds")
p1 <- plot_cell_trajectory(cds,
                           color_by = "Pseudotime", 
                           size = 1,
                           show_backbone = TRUE)
p2 <- plot_cell_trajectory(cds,
                           color_by = "State", 
                           size = 1,
                           show_backbone = TRUE)+
  scale_color_nejm()
p3 <- plot_cell_trajectory(cds,
                           color_by = "group",
                           size = 1,
                           show_backbone = TRUE)+
  scale_color_brewer(palette="Paired")
p4 <- plot_cell_trajectory(cds,
                           color_by = "sub_cell_type",
                           size = 1,
                           show_backbone = TRUE)+
  scale_color_npg()
p <- p1+p2+p3+p4
ggsave("1.pdf",p,height = 8,width = 8)

Time_genes <-top_n(deg, n = 1000, desc(p_val_adj)) %>% pull(gene) %>% as.character() 
Time_genes <- unique(Time_genes)
p <- plot_pseudotime_heatmap(
  cds[Time_genes,],
  num_clusters = 2,
  show_rownames = T,
  return_heatmap = T,
  hmcols =colorRampPalette(rev(brewer.pal(9, "PRGn")))(100)
)
ggsave("heatmap3.pdf",p,height = 5,width = 5)
clusters <- cutree(p$tree_row, k = 2) 
table(clusters)
genes <- names(clusters[clusters == 1])
write.table(genes,
            quote = F,
            row.names = F,
            col.names = F,
            file = "./gene/SMC_genes-1.txt")

# mapping umap----
library(monocle)
library(Seurat)
library(ggplot2)
pbmc <- readRDS("~/20240103_Atherosis/result/Fig_SMC/subset_seurat/SMC_anno.rds")
head(pbmc@reductions$umap@cell.embeddings)
HSMM <- readRDS("~/20240103_Atherosis/result/Fig_SMC/Supply/SMC_cds.rds")
head(HSMM@phenoData@data)
pbmc@meta.data$Pseudotime <- HSMM@phenoData@data$Pseudotime 
head(pbmc@meta.data)
p3 <- plot_cell_trajectory(HSMM, color_by = "Pseudotime")+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('blue','cyan','green','yellow','orange','red'))
p3
mydata<- FetchData(pbmc,vars = c("tSNE_1","tSNE_2","Pseudotime"))
p <- ggplot(mydata,aes(x = tSNE_1,y =tSNE_2,colour = Pseudotime))+
  geom_point_rast(size = 1)+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#324cd9',"#25cb78",'yellow','orange','red'))
p4 <- p + theme_bw() + theme(panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), 
                             axis.line = element_line(colour = "black"))
p4
outdir <- "~/20240103_Atherosis/result/Fig_SMC/Supply/映射tsne.pdf"
ggsave(filename = paste0(outdir,"umap_pesudotime2.pdf"),
       plot = p4,
       height = 4,
       width = 5)

# vlnplot----
library(msigdbr)
library(Seurat)
library(clusterProfiler)
library(ggplot2)
library(ggpubr)
setwd("~/20240103_Atherosis/v2/result/Fig_SMC/5-supply/12-score")

seurat_obj <- readRDS('~/20240103_Atherosis/result/Fig_SMC/subset_seurat/SMC_anno.rds')
C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)
a <- as.data.frame(unique(C5_gene_sets$gs_name))
selected_gene_sets <- C5_gene_sets %>% filter(gs_name %in% c( "GOBP_PHENOTYPIC_SWITCHING"))
seurat_obj <-
  AddModuleScore(seurat_obj, features = list(selected_gene_sets$gene_symbol), name = "PHENOTYPIC_SWITCHING")
p <- VlnPlot(
  seurat_obj,
  features = "PHENOTYPIC_SWITCHING1",
  group.by = "sub_cell_type",
  pt.size = 0,
  sort = T
) +
  NoLegend() +
  labs(title = "Genes Score") +
  stat_compare_means(method = "t.test")+
  theme(aspect.ratio = 0.3,
        axis.title.x = element_blank()) +
  scale_fill_manual(values = rep(c("#7FC97F","#BEAED4","#E64B35FF","#FDC086", "#FFFF99", "#386CB0",
                                   "#F0027F", "#BF5B17","#4dBBD5FF","#3C5488FF","#F39B7FFF","#91D1C2FF"),2))
p
ggsave("PHENOTYPIC_SWITCHING.pdf",p,height = 5,width = 3)
feature <- "WTAP"
seurat_obj <-
  AddModuleScore(seurat_obj, features = list(feature), name = "WTAP")
p <- VlnPlot(
  seurat_obj,
  features = "WTAP",
  group.by = "sub_cell_type",
  pt.size = 0,
  sort = T
) +
  NoLegend() +
  labs(title = "WTAP score") +
  theme(aspect.ratio = 0.3,
        axis.title.x = element_blank()) +
  stat_compare_means(method = "t.test")+
  scale_fill_manual(values = rep(c("#7FC97F","#BEAED4","#E64B35FF","#FDC086", "#FFFF99", "#386CB0",
                                   "#F0027F", "#BF5B17","#4dBBD5FF","#3C5488FF","#F39B7FFF","#91D1C2FF"),2))
ggsave("WTAP2.pdf",p,height = 5,width = 3)

# heatmap ------
library(ClusterGVis)
library(RColorBrewer)
library(viridis)
library(monocle) 
library(Seurat) 
library(dplyr) 
setwd("~/20240103_Atherosis/v2/result/Fig_SMC/5-supply/2-beautymonocle")

cds <- readRDS("~/20240103_Atherosis/result/Fig_SMC/Supply/SMC_cds.rds")
deg <- read_delim("~/20240103_Atherosis/result/Fig_SMC/Supply/SMC_monocle_deg.txt")
Time_genes <-top_n(deg, n = 1000, desc(p_val_adj)) %>% pull(gene) %>% as.character() #原来是ordergene，比较少，可以换成deg比较多
df <- plot_pseudotime_heatmap2(cds[Time_genes,],
                               num_clusters = 2,
                               show_rownames = F,
                               return_heatmap = F)
str(df)
PsPseudoGenes <- read_csv("~/20240103_Atherosis/v2/code/Fig_SMC/5-supply/10-PSintersectPseudoGenes/PsPseudoGenes.csv")
PsPseudoGenes <- PsPseudoGenes$x
C2gene <- read.table("~/20240103_Atherosis/result/Fig_SMC/Supply/gene/SMC_genes-2.txt")
PsPseudoGenes <- intersect(C2gene$V1,PsPseudoGenes)
m6Afeatures = c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
                "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
                "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
C5_gene_sets <- msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::select(gs_name, gene_symbol)
a <- as.data.frame(unique(C5_gene_sets$gs_name))
selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% c( "GOBP_ATP_METABOLIC_PROCESS"  )) 
ATPgenes <- unique(selected_gene_sets$gene_symbol)
C1gene <- read.table("~/20240103_Atherosis/result/Fig_SMC/Supply/gene/SMC_genes-1.txt")
ATPgenes <- intersect(ATPgenes,C1gene$V1)
marker <- c("ACTA2","MYH10","KLF4","KLF2") #"MYH11" ,"MYL9" ,"TAGLN"
genes <- c(PsPseudoGenes,m6Afeatures,ATPgenes,marker)
genes
pdf("1genes.pdf",height = 8,width = 7)
p <- plot_pseudotime_heatmap2(cds[genes,],
                              num_clusters = 2,
                              # cores = 1,
                              show_rownames = T,
                              return_heatmap = T,
                              # hmcols =colorRampPalette(rev(brewer.pal(9, "PRGn")))(100))
                              hmcols = colorRampPalette(c("#096ab3", "white", "#e2450c"))(100))
print(p)
dev.off()

gene <- c("ACTA2","MYH10","KLF4","KLF2","WTAP","FOS","JUN","JUNB","MYC",
          "KLF10","ATP5PF","COX5A","ATP5F1A","PGK1","PDGFA","CD46","NOTCH2","NOTCH3")# 挑一些基因展示
pdf(file = "4one-branch.pdf",height = 6,width = 7)
p3 <- visCluster(object = df,plot.type = "both",markGenes = gene,
                 pseudotime_col = c("#096ab3", "#e2450c"),
                 ht.col.list = list(col_range = seq(-2,2,length.out = 100),
                                    col_color = colorRampPalette(brewer.pal(9,"PRGn"))(100)))
print(p3)
dev.off()

# go----
setwd("~/20240103_Atherosis/result/Fig_SMC/Supply/go")
library(clusterProfiler)
library(org.Hs.eg.db) 
library(dplyr)

gene <- read.delim("~/20240103_Atherosis/result/Fig_SMC/Supply/gene/SMC_genes-2.txt", 
                   sep = ",", header = F)
head(gene)
bp <-
  enrichGO(
    gene$V1,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP", 
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term <- bp@result
write.table(
  term,
  file = "gene2_go.txt",
  row.names = F,
  quote = F,
  sep = "\t"
)
saveRDS(bp, "gene2_go.rds")
setwd("~/20240103_Atherosis/result/Fig_SMC/Supply/go")
library(ggplot2)
df <- read.delim("gene1_go.txt")
head(df)
pathways <- c(
  "aerobic respiration",
  "cellular respiration",
  "nucleoside triphosphate metabolic process",
  "ribonucleoside triphosphate metabolic process",
  "ATP metabolic process",
  "purine ribonucleoside triphosphate metabolic process",
  "purine nucleoside triphosphate metabolic process",
  "oxidative phosphorylation",
  "cytoplasmic translation",
  "respiratory electron transport chain",
  "protein folding",
  "tricarboxylic acid cycle",
  "glycolytic process"
)
term1 <- df[df$Description %in% pathways,]
term1$labelx=rep(0,nrow(term1))
term1$labely=seq(nrow(term1),1)
p <- ggplot(data = term1,
            aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#e15e5c",width = 0.6) +
  geom_text(aes(x=labelx, y=labely, label = term1$Description),size=3,hjust =0,color = "black")+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("p1.pdf", plot = p,
       height = 3.5, width = 5)
df <- read.delim("gene2_go.txt")
head(df)
pathways <- c(
  "smooth muscle cell proliferation",
  "positive regulation of smooth muscle cell proliferation",
  "vascular associated smooth muscle cell proliferation",
  "smooth muscle cell migration",
  "positive regulation of smooth muscle cell migration",
  "vascular associated smooth muscle cell migration",
  "smooth muscle cell differentiation",
  "extracellular matrix organization",
  "cell-matrix adhesion",
  "regulation of extracellular matrix assembly"
)
term1 <- df[df$Description %in% pathways,]
term1$labelx=rep(0,nrow(term1))
term1$labely=seq(nrow(term1),1)
p <- ggplot(data = term1,
            aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#3c86db",width = 0.6) +
  geom_text(aes(x=labelx, y=labely, label = term1$Description),size=3,hjust =0,color = "black")+
  theme_classic()+
  theme(axis.text.y = element_blank(),axis.line.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank(), axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave("p3.pdf", plot = p,
       height = 3.5, width = 5)

# pesudotime pathway----
library(monocle) 
library(Seurat) 
library(dplyr) 
library(ggsci)
setwd("~/20240103_Atherosis/v2/result/Fig_SMC/5-supply/9-pseudotimepathway")
outdir <- "~/20240103_Atherosis/v2/result/Fig_SMC/5-supply/9-pseudotimepathway"

seurat_obj <-
  readRDS("~/20240103_Atherosis/result/Fig_SMC/subset_seurat/SMC_anno.rds")
C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)
a <- as.data.frame(unique(C5_gene_sets$gs_name))
selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% c(
    "GOBP_POSITIVE_REGULATION_OF_SMOOTH_MUSCLE_CELL_MIGRATION",
    "GOBP_POSITIVE_REGULATION_OF_SMOOTH_MUSCLE_CELL_PROLIFERATION",
    "GOBP_POSITIVE_REGULATION_OF_VASCULAR_ASSOCIATED_SMOOTH_MUSCLE_CELL_MIGRATION",
    "GOBP_POSITIVE_REGULATION_OF_VASCULAR_ASSOCIATED_SMOOTH_MUSCLE_CELL_PROLIFERATION",
    "GOBP_SMOOTH_MUSCLE_CELL_MIGRATION",
    "GOBP_SMOOTH_MUSCLE_CELL_PROLIFERATION",
    "GOBP_SMOOTH_MUSCLE_CELL_MATRIX_ADHESION",
    "GOBP_REGULATION_OF_SMOOTH_MUSCLE_CELL_CHEMOTAXIS",
    "GOBP_CELL_MATRIX_ADHESION",
    "GOBP_PHENOTYPIC_SWITCHING",
    "GOBP_ATP_METABOLIC_PROCESS",
    "GOBP_OXIDATIVE_PHOSPHORYLATION",
    "GOBP_AEROBIC_RESPIRATION",
    "GOBP_RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",
    "GOBP_CELLULAR_RESPIRATION",
    "GOBP_POSITIVE_REGULATION_OF_GLYCOLYTIC_PROCESS",
    "GOBP_NUCLEOSIDE_TRIPHOSPHATE_METABOLIC_PROCESS"
  ))
selected_gene_sets
cells_rankings <- AUCell_buildRankings(as.matrix(seurat_obj@assays$RNA@data))
cells_rankings
geneSets<-lapply(unique(selected_gene_sets$gs_name),
                 function(x){selected_gene_sets$gene_symbol[selected_gene_sets$gs_name==x]})
View(geneSets)
names(geneSets) <- unique(selected_gene_sets$gs_name)
saveRDS(geneSets,"geneSets.rds")
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.05)
cells_AUC
saveRDS(cells_AUC,"auc.rds")
aucs <- getAUC(cells_AUC)
saveRDS(aucs,"aucs.RDS")
aucs <- t(aucs)
aucs <- as.data.frame(aucs)
cds <- readRDS("~/20240103_Atherosis/result/Fig_SMC/Supply/SMC_cds.rds")
cdsdata <- pData(cds)
b <- cbind(cdsdata,aucs)
pData(cds) <- b
name <- colnames(aucs)
for(i in name){
  p <- plot_cell_trajectory(cds,
                            color_by = i,
                            size = 1,
                            show_backbone = TRUE) + 
    scale_color_gradientn(colours = c("#2c0ab9", "#dd332f","#f8b072","#eef314"))
  ggsave(paste0(i,".pdf"),p,height = 6,width =6)
}


# cor-----
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library(magrittr)
library(ggplot2)
library(readr)
library(dplyr)
library(forcats)
library(stringr)
setwd("~/20240103_Atherosis/v2/result/Fig_SMC/5-supply/14-cor")
outdir <- "~/20240103_Atherosis/v2/result/Fig_SMC/5-supply/14-cor/"
dir.create(outdir, recursive = T)
selected_features <- c("ACTA2","COX5A","MYH10","ATP5PF","ATP5F1A","PGK1",
                       "KLF4","JUN","KLF2","MYC","KLF10","JUNB","FOS")
selected_genes <- "WTAP"
selected_cell_type <- "SMC"

adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/SMC/SMC/k_20/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
p <- filtered_adata %>%
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features) %>%
  mutate(
    feature_x = fct_relevel(feature_x, selected_genes),
    feature_y = fct_relevel(feature_y, selected_features)
  ) %>%
  catdotplot(
    x = feature_x,
    y = feature_y,
    size = -log10(p_value),
    color = estimate,
    title = selected_cell_type
  ) +
  coord_fixed() +
  scale_color_gradient2(
    low = "#2009af", mid = "white", high = "#af2934")+
  theme(
    panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.text = element_text(face = "italic")
  )
p
ggsave(
  file.path(outdir, "genes_cor1.pdf"),
  plot = p,
  height = 4,
  width = 4,
)

selected_genes <- c("WTAP")
adata <- read_csv("~/Atherosis_0723/20230204_Atherosis/result/SMC/SMC/k_20/correlation/GO:BP.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
selected_features <- c(     
  "GOBP_CELLULAR_RESPIRATION",
  "GOBP_AEROBIC_RESPIRATION",
  "GOBP_ATP_METABOLIC_PROCESS",
  "GOBP_OXIDATIVE_PHOSPHORYLATION",
  "GOBP_RIBONUCLEOSIDE_METABOLIC_PROCESS",
  "GOBP_ATP_BIOSYNTHETIC_PROCESS",
  "GOBP_PHENOTYPIC_SWITCHING"
)
p <- filtered_adata %>%
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features) %>%
  dplyr:: mutate(
    feature_x = fct_relevel(feature_x, selected_genes),
    feature_y = fct_relevel(feature_y, selected_features)
  ) %>%
  catdotplot(
    x = feature_x,
    y = feature_y,
    size = -log10(p_value),
    color = estimate,
    dot_scale = 7,
    title = selected_cell_type
  ) +
  scale_color_gradient2(
    low = "#2009af", mid = "white", high = "#af2934")+
  coord_fixed() +
  theme(
    panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.text.x = element_text(face = "italic")
  )
p  
ggsave(
  file.path(outdir, "term_cor1.pdf"),
  plot = p,
  height = 4,
  width = 4,
)

# pesudotime gene----
library(ClusterGVis)
library(RColorBrewer)
library(viridis)
library(monocle) 
library(Seurat)
library(dplyr) 
library(ggsci)
setwd("~/20240103_Atherosis/v2/result/Fig_SMC/5-supply/13-gene")

cds <- readRDS("~/20240103_Atherosis/result/Fig_SMC/Supply/SMC_cds.rds")
PsPseudoGenes <- read_csv("~/20240103_Atherosis/v2/code/Fig_SMC/5-supply/10-PSintersectPseudoGenes/PsPseudoGenes.csv")
PsPseudoGenes <- PsPseudoGenes$x
C2gene <- read.table("~/20240103_Atherosis/result/Fig_SMC/Supply/gene/SMC_genes-2.txt")
PsPseudoGenes <- intersect(C2gene$V1,PsPseudoGenes)
m6Afeatures = c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
                "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
                "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
C5_gene_sets <- msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::select(gs_name, gene_symbol)
a <- as.data.frame(unique(C5_gene_sets$gs_name))
selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% c( "GOBP_ATP_METABOLIC_PROCESS"  )) 
ATPgenes <- unique(selected_gene_sets$gene_symbol)
C1gene <- read.table("~/20240103_Atherosis/result/Fig_SMC/Supply/gene/SMC_genes-1.txt")
ATPgenes <- intersect(ATPgenes,C1gene$V1)
marker <- c("ACTA2","MYH10","KLF4","KLF2") 
genes <- c(PsPseudoGenes,m6Afeatures,ATPgenes,marker)
genes
C5_gene_sets <- msigdbr::msigdbr(species = "human", category = "C5") %>% dplyr::select(gs_name, gene_symbol)
a <- as.data.frame(unique(C5_gene_sets$gs_name))
selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% c( "GOBP_PHENOTYPIC_SWITCHING"  )) 
genes <- unique(selected_gene_sets$gene_symbol)
deg <- read_delim("~/20240103_Atherosis/result/Fig_SMC/Supply/SMC_monocle_deg.txt")
genes <- intersect(genes,deg$gene) # SOD2
s.genes <- c("WTAP","ACTA2","MYH10","JUN","FOS","JUNB","KLF4","ATP5PF","ENO1",
             "AGT","CCL2","KLF2","PPP1R15A","HIF1A","NDUFB2","PPP1CB","ATP5F1A",
             "EGR1","TSPO","SDHC","NDUFS4")
for(i in s.genes){
  p <- plot_genes_in_pseudotime(cds[i,], 
                                color_by = "sub_cell_type")+ 
    scale_color_npg()
  ggsave(paste0(i,".pdf"),p,height = 2,width = 4)
}

s.genes <- c("EGR1","JUN")
for(i in s.genes){
  p <- plot_genes_in_pseudotime(cds[i,], 
                                color_by = "sub_cell_type")+ 
    scale_color_npg()
  ggsave(paste0(i,".pdf"),p,height = 2,width = 4)
}

# pesudotime gene2 -----
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(AUCell)
library(SCENIC)
library(Seurat)
library(pheatmap)
library(patchwork)
library(qs)
library(stringr)
library(monocle)
library(dplyr)
library(tidyverse)
setwd("~/20240103_Atherosis/v2/result/Fig_SMC/5-supply/11-pesudotimegene/2")

seurat.obj <- readRDS("~/20240103_Atherosis/result/Fig_SMC/subset_seurat/SMC_anno.rds")
expr <- as.matrix(seurat.obj@assays$RNA@data)
expr <- as.data.frame(expr)
expr <- as.data.frame(t(expr))
expr[1:5,1:5]
cds <- readRDS("~/20240103_Atherosis/result/Fig_SMC/Supply/SMC_cds.rds")
cdsdata <- pData(cds)
expr2 <- subset(expr,rownames(expr) %in% rownames(cdsdata))
genes <- c("WTAP","MYH10","ACTA2","FOS","KLF2","JUNB","ATP5PF","NDUFB2","ENO1","CCL2","HIF1A","AGT")
expr3 <- expr2[,colnames(expr2) %in% genes]
expr3[1:5,1:5]
for( i in genes){
  b <- as.data.frame(expr3[,colnames(expr3) %in% i])
  colnames(b) <- paste0(i)
  df <- cbind(cdsdata,b)
  pData(cds) <- df
  p <-  plot_cell_trajectory(cds,
                             color_by = i,
                             size = 1,
                             show_backbone = TRUE) + 
    viridis::scale_color_viridis(option = "C")
  ggsave(paste0(i,".pdf"),p,height = 5,width = 5)
}
genes <- c("EGR1","JUN")
for( i in genes){
  b <- as.data.frame(expr3[,colnames(expr3) %in% i])
  colnames(b) <- paste0(i)
  df <- cbind(cdsdata,b)
  pData(cds) <- df
  p <-  plot_cell_trajectory(cds,
                             color_by = i,
                             size = 1,
                             show_backbone = TRUE) + 
    viridis::scale_color_viridis(option = "C")
  ggsave(paste0(i,".pdf"),p,height = 5,width = 5)
}



