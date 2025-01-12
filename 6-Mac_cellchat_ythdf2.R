# single --------
setwd("~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/YTHDF2/")
library(CellChat)
library(patchwork)
library(Seurat, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat)
library(tidyverse)
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/computing/custom_function.R")

seurat_obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
seurat_obj <- NormalizeData(seurat_obj)
selected_genes <- c("YTHDF2")
selected_cell_type <- "Macrophage"
outdir <- "~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/YTHDF2/"
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
seurat_obj[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj[[str_c(selected_genes, "_group")]]))

group <- "high"
output.dir <-
  paste0( "./",group, "/") 
dir.create(output.dir, recursive = T)
seurat.obj <- seurat_obj[, seurat_obj[["YTHDF2_group"]] == group] #
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
  paste0("./", group, "/") 
dir.create(output.dir, recursive = T)
seurat.obj <- seurat_obj[, seurat_obj[["YTHDF2_group"]] == group] 
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
saveRDS(cellchat, file = paste0(output.dir, "cellchat-low.rds"))

# multi --------
library(CellChat)
library(patchwork)
library(cowplot)
setwd("~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/YTHDF2/HighVSLow/")
outdir <- "~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/YTHDF2/HighVSLow/"

cellchat.control  <- readRDS("~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/YTHDF2/high/cellchat-high.rds")
cellchat.case  <- readRDS("~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/YTHDF2/low/cellchat-low.rds")
cellchat.control <- netAnalysis_computeCentrality(cellchat.control, slot.name = "netP")
cellchat.case <- netAnalysis_computeCentrality(cellchat.case, slot.name = "netP")
object.list <- list(High = cellchat.control, 
                    Low = cellchat.case)
cellchat <- mergeCellChat(
  object.list,
  add.names = names(object.list))
cellchat
gg1 <- netVisual_heatmap(cellchat,font.size = 5,font.size.title = 5,color.heatmap = c("#7A577A", "#ea9241"))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",font.size = 5,font.size.title = 5,color.heatmap = c("#7A577A", "#ea9241"))

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) 
gg <- list()
pdf("6-2D.pdf",height = 5,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("3-热图.pdf",width = 8,height = 4)
p <- gg1 + gg2
print(p)
dev.off()
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
dev.off()
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = c("GRN","TNF","HSPG","TIGIT","NECTIN","THY1","IGF","CD34","PERIOSTIN","CD6","ALCAM","CDH","VCAM","CXCL","ITGB2"),
                                        title = names(object.list)[i],width = 5, height = 7.5,font.size = 4,font.size.title = 6, color.heatmap = "YlGnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = c("GRN","TNF","HSPG","TIGIT","NECTIN","THY1","IGF","CD34","PERIOSTIN","CD6","ALCAM","CDH","VCAM","CXCL","ITGB2"),
                                        title = names(object.list)[i+1], width = 5, height = 7.5,font.size = 4,font.size.title = 6, color.heatmap = "YlGnBu")
pdf("7-heatmap3-3.pdf",height = 4,width = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
levels(object.list[[1]]@idents)
table(object.list$High@meta)
netVisual_bubble(cellchat,
                 sources.use = c(3),  
                 targets.use = c(1:2,4:15),  
                 font.size = 6,font.size.title = 6,
                 remove.isolate = T,
                 signaling = c("CXCL","TNF","VCAM","NECTIN","ALCAM","ITGB2","CD6","GRN"),
                 comparison = c(1:2), angle.x = 45)  
ggsave("6-气泡图1.pdf",height = 4,width = 6)
netVisual_bubble(cellchat,
                 sources.use = c(1:2,4:15),  
                 targets.use = c(3),  
                 remove.isolate = T,
                 font.size = 6,font.size.title = 6,
                 signaling = c("CXCL","TNF","VCAM","NECTIN","ALCAM","ITGB2","CD6","GRN","ADGRE5"),
                 comparison = c(1:2), angle.x = 45)  
ggsave("6-气泡图3.pdf",height = 3,width = 6)
pathways.show <- c("TNF")
weight.max <- getMaxWeight(object.list[1],slot.name = c('netP') ,attribute =pathways.show)
pdf("TNF.pdf",height = 5,width = 5)
netVisual_aggregate(object.list[[1]],signaling = pathways.show,layout = "chord",
                    edge.weight.max = weight.max[1],edge.width.max = 10,
                    signaling.name = paste(pathways.show,names(object.list)[1]),
                    pt.title = 4,title.space = 4,
)
dev.off()
pdf("9-通路基因_TNF.pdf",height = 8,width = 8)
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,3:15), targets.use = c(1,2,3:15),  slot.name = "net",lab.cex = 1,
                     title.name = paste0("TNF signaling from all celltypes - ", names(object.list)[1]),
                     signaling = c("TNF"),big.gap = 5,small.gap = 1,annotationTrackHeight = c(0.05))
dev.off()
pdf("9-通路基因_VCAM.pdf",height = 8,width = 8)
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,3:15), targets.use = c(1,2,3:15),  slot.name = "net",lab.cex = 1,
                     title.name = paste0("VCAM signaling from all celltypes - ", names(object.list)[1]),
                     signaling = c("VCAM"),big.gap = 5,small.gap = 1,annotationTrackHeight = c(0.05))
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
pdf("10-通路基因_CXCL.pdf",height = 8,width = 8)
netVisual_chord_gene(object.list[[1]], sources.use = c(1,2,3:15), targets.use = c(1,2,3:15), 
                     signaling = c("CXCL"),slot.name = 'net', net = net.down, lab.cex = 1, small.gap = 3.5, 
                     title.name = paste0("Upregualted CXCL signaling in YTHDF2-", names(object.list)[1])) 
dev.off()
save(object.list, file = "/home/pingxr/20240103_Atherosis/v3/result/Fig_Mac/cellchat/YTHDF2/cellchat_object.list.RData")
save(cellchat, file = "/home/pingxr/20240103_Atherosis/v3/result/Fig_Mac/cellchat/YTHDF2/cellchat_merged_.RData")

#gene expression/cor-------
setwd("~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/YTHDF2/ythdf2-geneexpr/")
outdir <- "~/20240103_Atherosis/v3/result/Fig_Mac/cellchat/YTHDF2/ythdf2-geneexpr/"
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
selected_genes <- c("YTHDF2")
selected_cell_type <- "Macrophage"
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
seurat_obj[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj[[str_c(selected_genes, "_group")]]))
genes <- c("CXCL8","ITGB2","TNF","TNFRSF1B","ITGB1","ITGA4")
df <- data_frame()
for(i in genes){
  p <- VlnPlot(
    seurat_obj,
    features = i,
    group.by = "YTHDF2_group",
    assay = "RNA"
  )
  data <- p$data
  data$gene <- paste0(i)
  colnames(data) <- c("value","group","gene")
  df <- rbind(df,data)
}
library(ggridges)
cols=c("#ea9241","#7A577A", "#EE757A","#A78CC1", "#DBCA85","#FFAE5D","#78ADD5","#CE9D7D","#84C780","#DED2E6","#C8E2EF")
ggplot(df, aes(x = value, y = gene,fill=group,color=group)) +
  geom_density_ridges(
    alpha = 0.6,          
    color = 'white',      
    rel_min_height = 0.01, 
    scale = 1,          
    quantile_lines = TRUE, 
    quantiles = 2,       ）
    size = 0.5,           
    show.legend = TRUE,  
    bandwidth = 0.2
  ) +
  scale_fill_manual(values = cols) + 
  scale_color_manual(values = cols)+ 
  theme_bw() +
  ylab("Gene")+
  xlab("Expression")+
  scale_x_continuous(limits = c(0.5, 5), breaks = seq(0.5, 5, by = 1)) +
  theme(axis.title.x = element_text(size = 8,color = "black"),  
        axis.title.y = element_text(size = 8,color = "black"), 
        axis.text.x = element_text(size = 8,color = "black"),   
        axis.text.y = element_text(size = 8,color = "black"))+
  theme(plot.title = element_text(hjust = 0.2, vjust = 2,size = 5)  
  )
ggsave("all.pdf",height = 3,width = 4.5)

source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
selected_genes <- "YTHDF2"
selected_cell_type <- "Macrophage"
genes <- c(  "ITGA4","TNFRSF1B","CXCL2","CD6", 
             "ITGB1","CXCR4"  ,"CXCL8","ITGB2","TNF",
             "TNFRSF1A", 
             "ICAM1","ICAM2",
             "CXCL3","CXCL12"
)
selected_features <- unique(genes)
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

#pathway score-------
setwd("~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/4-ythdf2-termcor")
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
x <- "YTHDF2"
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
  "GOBP-CYTOKINE-PRODUCTION",
  "GOBP-CYTOKINE-PRODUCTION-INVOLVED-IN-INFLAMMATORY-RESPONSE",
  "GOBP-CELL-CELL-ADHESION",
  "GOBP-REGULATION-OF-CELL-ADHESION",
  "GOBP-LEUKOCYTE-CELL-CELL-ADHESION",
  "GOBP-INFLAMMATORY-RESPONSE"
)

df <- data_frame()
for(i in terms){
  p <- VlnPlot(
    sub_seurat_obj,
    features = paste0(i),
    group.by = "YTHDF2_group",
    assay = "RNA"
  )
  p
  data <- p$data
  data$term <- paste0(i)
  colnames(data) <- c("value","group","term")
  df <- rbind(df,data)
}
table(df$term)
table(df$group)
df$term <- gsub("GOBP-","",df$term)
df$term <- gsub("-"," ",df$term)
df$term <- tolower(df$term)
p <- ggplot(df,aes(x = term,y = value,fill = group)) +
  geom_split_violin(alpha = 0.8, trim = F,color = NA,width = 1) +
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,color = 'black',hjust = 1,size = 5),
        legend.position = 'top') +
  scale_fill_manual(values = c("#eaa86c","#ab7cab"))+
  stat_compare_means(aes(group=group),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif",
                     label.y = max(df$value), size = 3)
ggsave("all.pdf",p,height = 5,width = 6)





