# single-----
setwd("~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/METTL3/")
library(CellChat)
library(patchwork)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat)
library(tidyverse)
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/computing/custom_function.R")

seurat_obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
seurat_obj <- NormalizeData(seurat_obj)
selected_genes <- c("METTL3")
selected_cell_type <- "Macrophage"
outdir <- "~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/METTL3/"
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
seurat_obj[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj[[str_c(selected_genes, "_group")]]))
group <- "high"
output.dir <-
  paste0( "./",group, "/") # must have "/"
dir.create(output.dir, recursive = T)
seurat.obj <- seurat_obj[, seurat_obj[["METTL3_group"]] == group] #
names(table(Idents(seurat.obj)))
print(table(Idents(seurat.obj)))
expr <- seurat.obj@assays$RNA@data
data.input <- expr
dim(data.input)
data.input[1:4, 1:4]
meta <- as.data.frame(Idents(seurat.obj))
colnames(meta) <- "labels"
unique(meta$labels) 
cellchat <-
  createCellChat(object = data.input,
                 meta = meta,
                 group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <-
  setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) 
groupSize <-
  as.numeric(table(cellchat@idents)) 
CellChatDB <-
  CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <-
  subsetData(cellchat) 
future::plan("multisession", workers = 10) 
cellchat <-
  identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
saveRDS(cellchat, file = paste0(output.dir, "cellchat-high.rds"))

group <- "low"
output.dir <-
  paste0("./",group, "/") # must have "/"
dir.create(output.dir, recursive = T)
seurat.obj <- seurat_obj[, seurat_obj[["METTL3_group"]] == group] #
names(table(Idents(seurat.obj)))
print(table(Idents(seurat.obj)))
expr <- seurat.obj@assays$RNA@data
data.input <- expr
dim(data.input)
data.input[1:4, 1:4]
meta <- as.data.frame(Idents(seurat.obj))
colnames(meta) <- "labels"
unique(meta$labels) # check the cell labels
cellchat <-
  createCellChat(object = data.input,
                 meta = meta,
                 group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <-
  setIdent(cellchat, ident.use = "labels") 
levels(cellchat@idents) 
groupSize <-
  as.numeric(table(cellchat@idents)) 
CellChatDB <-
  CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <-
  subsetData(cellchat) 
future::plan("multisession", workers = 10) 
cellchat <-
  identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
saveRDS(cellchat, file = paste0(output.dir, "cellchat-low.rds"))

# multi----
library(CellChat)
library(patchwork)
library(cowplot)
setwd("~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/METTL3/HighVSLow/")
outdir <- "~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/METTL3/HighVSLow/"

cellchat.control  <- readRDS("~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/METTL3/high/cellchat-high.rds")
cellchat.case  <- readRDS("~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/METTL3/low/cellchat-low.rds")
cellchat.control <- netAnalysis_computeCentrality(cellchat.control, slot.name = "netP")
cellchat.case <- netAnalysis_computeCentrality(cellchat.case, slot.name = "netP")
object.list <- list(High = cellchat.control, 
                    Low = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list))
cellchat
library(CellChat)
cellchat <- updateCellChat(cellchat)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf("4-网络图.pdf",height = 5,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count,
                   weight.scale = T,
                   label.edge= F,
                   edge.weight.max = weight.max[2],
                   edge.width.max = 5, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

library(ComplexHeatmap)
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways,object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = c("TNF","CD40","VCAM","SELPLG","PECAM1","EGF","GRN","BAFF","ADGRE5"), 
                                        title = names(object.list)[i],width = 5, height = 7.5,font.size = 4,font.size.title = 6, color.heatmap = "YlGn")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = c("TNF","CD40","VCAM","SELPLG","PECAM1","EGF","GRN","BAFF","ADGRE5"), 
                                        title = names(object.list)[i+1], width = 5, height = 7.5,font.size = 4,font.size.title = 6, color.heatmap = "YlGn")
pdf("6-heatmap5.pdf",height = 4,width = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

table(object.list$High@meta)
netVisual_bubble(cellchat,
                 sources.use = c(1:15),  
                 targets.use = c(3),
                 font.size = 6,font.size.title = 6,
                 comparison = c(1:2), angle.x = 45,remove.isolate = T,
                 signaling = c("TNF","PECAM1","VCAM","SELPLG","EGF","CD40","GRN","BAFF","ADGRE5"))   
ggsave("8-气泡图6.pdf",height = 4,width = 6.5)
netVisual_bubble(cellchat,
                 sources.use = c(3),  
                 targets.use = c(1:15),
                 font.size = 6,font.size.title = 6,
                 comparison = c(1:2), angle.x = 45,remove.isolate = T,
                 signaling = c("TNF","PECAM1","VCAM","SELPLG","EGF","CD40","GRN","BAFF","ADGRE5"))   
ggsave("8-气泡图7.pdf",height = 4,width = 6.5)

pathways.show <- c("CD45")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("CD45.pdf",height = 15,width = 10)
par(mfrow = c(1,2), xpd=TRUE)  
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()
pathways.show <- c("CD40")
weight.max <- getMaxWeight(object.list[1],slot.name = c('netP') ,attribute =pathways.show)
pdf("CD40.pdf",height = 5,width = 5)
netVisual_aggregate(object.list[[1]],signaling = pathways.show,layout = "chord",
                    edge.weight.max = weight.max[1],edge.width.max = 10,
                    signaling.name = paste(pathways.show,names(object.list)[1]),
                    pt.title = 4,title.space = 4,
)
dev.off()

pdf("9-通路基因_CD45.pdf",height = 15,width = 15)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3:15), targets.use = c(1,2,3:15),  slot.name = "net",lab.cex = 1,
                       title.name = paste0("CD45 signaling from all celltypes - ", names(object.list)[i]),
                       signaling = c("CD45"),big.gap = 5,small.gap = 1,annotationTrackHeight = c(0.05))
}
dev.off()
pdf("9-通路基因_GRN.pdf",height = 8,width = 8)
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,3:15), targets.use = c(1,2,3:15),  slot.name = "net",lab.cex = 1,
                     title.name = paste0("GRN signaling from all celltypes - ", names(object.list)[1]),
                     signaling = c("GRN"),big.gap = 5,small.gap = 1,annotationTrackHeight = c(0.05))
dev.off()

pos.dataset = "Low"  
features.name = pos.dataset
cellchat <-
  identifyOverExpressedGenes(
    cellchat,
    group.dataset = "datasets",
    pos.dataset = pos.dataset,
    features.name = features.name,
    only.pos = FALSE,
    thresh.pc = 0.1,
    thresh.fc = 0.1,
    thresh.p = 1
  )
net <- netMappingDEG(cellchat, features.name = features.name)
net.up <- subsetCommunication(cellchat, net = net, datasets = "Low",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "High",ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pdf("10-通路基因_PECAM1.pdf",height = 8,width = 8)
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,3:15), targets.use = c(1,2,3:15), 
                     signaling = c("PECAM1"),slot.name = 'net', net = net.down, lab.cex = 1, small.gap = 3.5, 
                     title.name = paste0("Upregualted PECAM1 signaling in METTL3-", names(object.list)[1])) 
dev.off()
save(object.list, file = "~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/METTL3/cellchat_object.list.RData")
save(cellchat, file = "~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/METTL3/cellchat_merged_.RData")


# gene expression/ cor----
setwd("~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/METTL3/mettl3-geneexpr/")
outdir <- "~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/METTL3/mettl3-geneexpr/"
library(CellChat)
library(patchwork)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat)
library(tidyverse)
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/computing/custom_function.R")

seurat_obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- subset(seurat_obj,cell_type == "Macrophage")
selected_genes <- c("METTL3")
selected_cell_type <- "Macrophage"
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
seurat_obj[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj[[str_c(selected_genes, "_group")]]))
genes <- c("CD40LG","ITGA2B","ITGB3","ITGA5","ITGB1","ITGAM","ITGB2","CD40", #CD40
           "TNF","TNFRSF1A","TNFRSF1B","TNFSF13B","TNFRSF17","TNFRSF13C", #TNF
           "VCAM1","ITGA4","ITGB1", #VCAM
           "SELPLG","SELL","SELP", #SELPLG
           "GRN","SORT1",#GRN
           "PECAM1" #PECAM1
)
library(gghalves)
for(i in genes){
  p <- VlnPlot(
    seurat_obj,
    features = i,
    group.by = "METTL3_group",
    assay = "RNA"
  )
  data <- p$data
  data$gene <- paste0(i)
  colnames(data) <- c("value","group","gene")
  
  p <-ggplot(data = data) +
    geom_half_violin(
      aes(x = gene, y = value, fill = group,color = group),
      data = data %>% filter(group == "high"),
       side = "l"
    ) +
    geom_half_violin(
      aes(x = gene, y = value, fill = group,color = group),
      data = data %>% filter(group == "low"),
       side = "r"
    ) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_blank()) +
    xlab("") +
    ylab("expression") +
    geom_point(aes(x = gene, y = value, fill = group),
               stat = 'summary', fun = mean,
               position = position_dodge(width = 0.2)) +
    stat_summary(aes(x = gene, y = value, fill = group),
                 fun.min = function(x) {quantile(x)[2]},
                 fun.max = function(x) {quantile(x)[4]},
                 geom = 'errorbar', color = 'black',
                 width = 0.01, size = 0.5,
                 position = position_dodge(width = 0.2)) +
    stat_compare_means(
      aes(x = gene, y = value, fill = group),
      data = data,
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                         symbols = c("***", "**", "*", "")),
      label = "p.signif",
      method = "t.test"
    ) +
    theme(axis.text.x = element_text(angle = 0, hjust = 1),  
          legend.position = "top") +
    scale_fill_manual(values = c("low" = '#4c659f', "high" = '#fe6958'))+
    scale_color_manual(values = c("low" = '#4c659f', "high" = '#fe6958'))
  ggsave(paste0(i,".pdf"),height = 3,width = 2)
}

source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
selected_genes <- "METTL3"
selected_cell_type <- "Macrophage"
selected_features <- c("CD40LG","ITGA2B","ITGB3","ITGA5","ITGB1","ITGAM","ITGB2",#"CD40", #CD40
                       "TNF","TNFRSF1A","TNFRSF1B","TNFSF13B","TNFRSF17","TNFRSF13C", #TNF
                       "VCAM1","ITGA4","ITGB1", #VCAM
                       "SELPLG","SELL","SELP", #SELPLG
                       "GRN","SORT1",#GRN
                       "PECAM1" #PECAM1
)
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
  theme(
    panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.text = element_text(face = "italic")
  )
p   
ggsave("genes_cor_dotplot.pdf",
       plot = p,
       height = 3,
       width = 3,
)

# pathway score-----
setwd("~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/4-mettl3-termcor")
library(Seurat)
library(tidyverse) 
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library(introdataviz)

adata <-readRDS( "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/GO:BP.scores.rds" )
seurat_obj <-readRDS( "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/data/seurat_obj.rds" )
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "Macrophage")
pacman::p_load(Seurat)
pacman::p_load(tidyverse)
adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj
x <- "METTL3"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c(x, "_group")]]))

sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
a <- as.data.frame(rownames(sub_seurat_obj))
terms <- c(
  "GOBP-CYTOKINE-MEDIATED-SIGNALING-PATHWAY",
  "GOBP-MACROPHAGE-CYTOKINE-PRODUCTION",
  "GOBP-CYTOKINE-PRODUCTION",
  "GOBP-CELL-CELL-ADHESION",
  "GOBP-REGULATION-OF-CELL-ADHESION",
  "GOBP-LEUKOCYTE-CELL-CELL-ADHESION",
  "GOBP-REGULATION-OF-ANTIGEN-RECEPTOR-MEDIATED-SIGNALING-PATHWAY",
  "GOBP-RESPONSE-TO-LECTIN",
  "GOBP-CELL-CELL-SIGNALING",
  "GOBP-ERBB-SIGNALING-PATHWAY",
  "GOBP-ERBB2-SIGNALING-PATHWAY",
  "GOBP-REGULATION-OF-ERBB-SIGNALING-PATHWAY",
  "GOBP-CELL-CYCLE"
)
for(i in terms){
  p <- VlnPlot(
    sub_seurat_obj,
    features = paste0(i),
    group.by = "METTL3_group",
    assay = "RNA"
  )
  p
  data2 <- p$data
  a <- tolower(paste0(i))
  a <- gsub("gobp-","",a)
  a <- gsub("-"," ",a)
  data2$term <- a
  colnames(data2) <- c("value","group","term")
  p <- ggplot(data2,aes(x = term,y = value,fill = group)) +
    geom_split_violin(alpha = 0.6, trim = F,color = NA,width = 1) +
    stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
    stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
                 size = 0.3,
                 position = position_dodge(0.2)) +
    theme_bw(base_size = 5) +
    theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1),
          legend.position = 'top') +
    scale_fill_manual(values = c("#fa6a56", "#0570b1"))+
    stat_compare_means(aes(group=group),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "NS")),label = "p.signif",
                       label.y = max(data2$value), size = 3)
  
  ggsave(paste0(a,".pdf"),height = 3,width = 1.5)
}





