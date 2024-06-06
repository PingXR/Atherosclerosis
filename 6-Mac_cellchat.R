# mettl3 cellchat----
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
  paste0( "./",group, "/")
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
pdf(file = paste0(output.dir, "netVisual_circle_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
pdf(file = paste0(output.dir, "netVisual_circle_by_weight.pdf"))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
)
dev.off()
pdf(file = paste0(output.dir, "netVisual_circle_by_cell_type.pdf"),width = 10,height = 10)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <-
    matrix(
      0,
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
  mat2[i,] <- mat[i,]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = T,
    arrow.width = 0.2,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()
pdf(file = paste0(output.dir, "netVisual_chord_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
group <- "low"
output.dir <-
  paste0("./",group, "/")
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
saveRDS(cellchat, file = paste0(output.dir, "cellchat-low.rds"))
pdf(file = paste0(output.dir, "netVisual_circle_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
pdf(file = paste0(output.dir, "netVisual_circle_by_weight.pdf"))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
)
dev.off()
pdf(file = paste0(output.dir, "netVisual_circle_by_cell_type.pdf"),width = 10,height = 10)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <-
    matrix(
      0,
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
  mat2[i,] <- mat[i,]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = T,
    arrow.width = 0.2,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()
pdf(file = paste0(output.dir, "netVisual_chord_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
saveRDS(cellchat, file = paste0(output.dir, "cellchat-low.rds"))
# mettl3 mutilcellchat-----
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
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
pdf("5-2D.pdf",height = 5,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
dev.off()
library(ComplexHeatmap)
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways,object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i],width = 5, height = 7.5,font.size = 4,font.size.title = 6, color.heatmap = "YlGn")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 7.5,font.size = 4,font.size.title = 6, color.heatmap = "YlGn")
pdf("6-heatmap3.pdf",height = 10,width = 7)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
table(object.list$High@meta)
netVisual_bubble(cellchat,
                 sources.use = c(1:6,9,11,14,15),  
                 targets.use = c(1:6,9,11,14,15),
                 font.size = 6,font.size.title = 6,
                 comparison = c(1:2), angle.x = 45,remove.isolate = T,
                 signaling = c("TNF","PECAM1","VCAM","SELPLG","CD45","EGF"))   #comparison = c(1, 2)就是group里面有几组 2组 Before Late,
ggsave("8-气泡图4.pdf",height = 4,width = 6.5)
pathways.show <- c("EGF")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("EGF.pdf",height = 15,width = 10)
par(mfrow = c(1,2), xpd=TRUE) 
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()
pathways.show <- c("SELPLG")
weight.max <- getMaxWeight(object.list[1],slot.name = c('netP') ,attribute =pathways.show)
pdf("SELPLG.pdf",height = 5,width = 5)
netVisual_aggregate(object.list[[1]],signaling = pathways.show,layout = "chord",
                    edge.weight.max = weight.max[1],edge.width.max = 10,
                    signaling.name = paste(pathways.show,names(object.list)[1]),
                    pt.title = 4,title.space = 4,
)
dev.off()
save(object.list, file = "~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/METTL3/cellchat_object.list.RData")
save(cellchat, file = "~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/METTL3/cellchat_merged_.RData")
# mettl3 genecor----
seurat_obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- subset(seurat_obj,cell_type == "Macrophage")
selected_genes <- c("METTL3")
selected_cell_type <- "Macrophage"
outdir <- "~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/3-mettl3-geneexpr/"
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
seurat_obj[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj[[str_c(selected_genes, "_group")]]))
genes <- c("SELL","SELP","SELE","EGFR","ERBB2","AREG","ITGB7","VCAM1","ITGA9")
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
  p <- ggplot()+
    geom_half_violin(
      data = data %>% filter(group == "high"),
      aes(x = gene,y = value),colour="#e64b35ff",fill="#e64b35ff",side = "l"
    )+
    geom_half_violin(
      data = data %>% filter(group == "low"),
      aes(x = gene,y = value),colour='#00AFBB',fill='#00AFBB',side = "r"
    )+
    theme_bw()+
    theme(panel.grid=element_blank(),axis.text.x=element_blank())+
    xlab("")+
    ylab("expression")+
    geom_point(data = data, aes(x = gene,y = value, fill = group),
               stat = 'summary', fun=mean,
               position = position_dodge(width = 0.2))+
    stat_summary(data = data, aes(x = gene,y = value, fill = group),
                 fun.min = function(x){quantile(x)[2]},
                 fun.max = function(x){quantile(x)[4]},
                 geom = 'errorbar', color='black',
                 width=0.01,size=0.5,
                 position = position_dodge(width = 0.2) )+
    stat_compare_means(data = data, aes(x = gene,y = value),
                       symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "-")),
                       label = "p.signif",
                       label.y = max(data$value),
                       paired = "wilcox.test", 
                       hide.ns = F)+
    theme(axis.text.x = element_text(angle = 0, hjust = 1),  
          legend.position = "top",legend.key = element_rect(fill = c("#00AFBB","#e64b35ff")),
          legend.justification = "right")
  ggsave(paste0(i,".pdf"),height = 3,width = 2)
}
source("custom_plot_function.R")
adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
selected_genes <- "METTL3"
selected_cell_type <- "Macrophage"
selected_features <- c("ITGA4","ITGB1","PECAM1","SELPLG","TNF","PTPRC","TNFRSF1A","TNFRSF1B")
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
# mettl3 termcor----
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

# --------#m6A-相关性------
selected_genes <- c("METTL3") 
# Ready to analyze cell type?
selected_cell_type <- "Macrophage"
#### 气泡图展示关注基因与通路的相关性 
adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/GO:BP.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
selected_features <- c(  
  "GOBP_RESPONSE_TO_LECTIN",
  "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
  "GOBP_REGULATION_OF_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_CELL_CYCLE",
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_CYTOKINE_PRODUCTION",
  "GOBP_REGULATION_OF_CELL_ADHESION",
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_CELL_CELL_SIGNALING",
  "GOBP_ERBB_SIGNALING_PATHWAY",
  "GOBP_ERBB2_SIGNALING_PATHWAY",
  "GOBP_REGULATION_OF_ERBB_SIGNALING_PATHWAY",
  "GOBP_INFLAMMATORY_RESPONSE"
)
df <- filtered_adata %>%
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features) 
df$feature_y <- gsub("GOBP_", "", df$feature_y)
df$feature_y <- gsub("_", " ", df$feature_y)
df$feature_y <- tolower(df$feature_y)
p <- df %>%
  dplyr:: mutate(
    feature_x = fct_relevel(feature_x, selected_genes),
    feature_y = fct_relevel(feature_y, selected_features)
  ) %>%
  ggplot(
    aes(x = feature_x,
        y = feature_y)
  ) +
  geom_point(aes(size=-log10(p_value), color=estimate),
             shape = 19 
             # stroke = 3
  ) +
  scale_color_gradient2(low = "#0045d1", high = "red", mid = "yellow") +
  coord_fixed() +
  theme_bw()+
  theme(
    panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.text.x = element_text(face = "italic"),
    axis.line = element_line(color = "black",size = 0.2)
  ) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5)) 
ggsave(
  "cor.pdf",
  plot = p,
  height = 6,
  width = 8,
)






# ythdf2 cellchat----
library(CellChat)
library(patchwork)
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
pdf(file = paste0(output.dir, "netVisual_circle_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
pdf(file = paste0(output.dir, "netVisual_circle_by_weight.pdf"))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
)
dev.off()
pdf(file = paste0(output.dir, "netVisual_circle_by_cell_type.pdf"),width = 10,height = 10)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <-
    matrix(
      0,
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
  mat2[i,] <- mat[i,]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = T,
    arrow.width = 0.2,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()
pdf(file = paste0(output.dir, "netVisual_chord_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
group <- "low"
output.dir <-
  paste0("./", group, "/") 
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
saveRDS(cellchat, file = paste0(output.dir, "cellchat-low.rds"))
pdf(file = paste0(output.dir, "netVisual_circle_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
pdf(file = paste0(output.dir, "netVisual_circle_by_weight.pdf"))
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Interaction weights/strength"
)
dev.off()
pdf(file = paste0(output.dir, "netVisual_circle_by_cell_type.pdf"),width = 10,height = 10)
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <-
    matrix(
      0,
      nrow = nrow(mat),
      ncol = ncol(mat),
      dimnames = dimnames(mat)
    )
  mat2[i,] <- mat[i,]
  netVisual_circle(
    mat2,
    vertex.weight = groupSize,
    weight.scale = T,
    arrow.width = 0.2,
    edge.weight.max = max(mat),
    title.name = rownames(mat)[i]
  )
}
dev.off()
pdf(file = paste0(output.dir, "netVisual_chord_by_count.pdf"))
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = T,
  label.edge = F,
  title.name = "Number of interactions"
)
dev.off()
saveRDS(cellchat, file = paste0(output.dir, "cellchat-low.rds"))

# ythdf2 mutilcellchat----
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
saveRDS(cellchat, file = "./merged_cellchat.rds")
gg1 <- netVisual_heatmap(cellchat,font.size = 5,font.size.title = 5,color.heatmap = c("#7A577A", "#ea9241"))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",font.size = 5,font.size.title = 5,color.heatmap = c("#7A577A", "#ea9241"))
pdf("3-热图.pdf",width = 8,height = 4)
p <- gg1 + gg2
print(p)
dev.off()
gg1 <- rankNet(cellchat,
               mode = "comparison", stacked = T,
               do.stat = TRUE) +coord_flip() +scale_fill_manual(values = c("#ab7cab", "#eaa86c"))
gg2 <- rankNet(cellchat,
               mode = "comparison", stacked = F, 
               do.stat = TRUE) +coord_flip()+scale_fill_manual(values = c("#ab7cab", "#eaa86c"))
gg1+gg2
ggsave("5-堆叠柱形图.pdf", plot = gg1+gg2, width = 12,height = 10)
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
library(ComplexHeatmap)
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways,object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i],width = 5, height = 7.5,font.size = 4,font.size.title = 6, color.heatmap = "YlGnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 7.5,font.size = 4,font.size.title = 6, color.heatmap = "YlGnBu")
pdf("7-heatmap3-2.pdf",height = 10,width = 7)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()
table(object.list$High@meta)
netVisual_bubble(cellchat,
                 sources.use = c(3), 
                 targets.use = c(1:2,4:15), 
                 font.size = 6,font.size.title = 6,
                 remove.isolate = T,
                 signaling = c("CXCL","TNF","VCAM","NECTIN","ALCAM","ITGB2","CD6","GRN"),
                 comparison = c(1:2), angle.x = 45)   
ggsave("6-气泡图1.pdf",height = 4,width = 6)
pathways.show <- c("ITGB2")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("ITGB2.pdf",height = 15,width = 10)
par(mfrow = c(1,2), xpd=TRUE)  
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]))
}
dev.off()
pathways.show <- c("CD6")
weight.max <- getMaxWeight(object.list[1],slot.name = c('netP') ,attribute =pathways.show)
pdf("CD6.pdf",height = 5,width = 5)
netVisual_aggregate(object.list[[1]],signaling = pathways.show,layout = "chord",
                    edge.weight.max = weight.max[1],edge.width.max = 10,
                    signaling.name = paste(pathways.show,names(object.list)[1]),
                    pt.title = 4,title.space = 4,
)
dev.off()
save(object.list, file = "~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/cellchat2/YTHDF2/cellchat_object.list.RData")
# ytdhf2 genecor----
seurat_obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- subset(seurat_obj,cell_type == "Macrophage")
selected_genes <- c("YTHDF2")
selected_cell_type <- "Macrophage"
outdir <- "~/20240103_Atherosis/v2/result/Fig_Mac/5-supply/3-ythdf2-geneexpr/"
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
seurat_obj[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj[[str_c(selected_genes, "_group")]]))
library(ggridges)
cols=c("#EE757A","#A78CC1", "#DBCA85","#FFAE5D","#78ADD5","#CE9D7D","#84C780","#DED2E6","#C8E2EF")
ggplot(df, aes(x = value, y = gene,fill=group,color=group)) + 
  geom_density_ridges(
    alpha = 0.6,          
    color = 'white',      
    rel_min_height = 0.01, 
    scale = 1,         
    quantile_lines = TRUE, 
    quantiles = 2,        
    size = 0.5,             
    show.legend = TRUE,  
    bandwidth = 0.2
  ) +
  scale_fill_manual(values = cols) + 
  scale_color_manual(values = cols)+
  theme_bw() +
  ylab("Gene")+
  xlab("Expression")+
  scale_x_continuous(limits = c(0.5, 5), breaks = seq(0.5, 5, by = 1))
ggsave("all2.pdf",height = 3,width = 5)
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
selected_genes <- "YTHDF2"
selected_cell_type <- "Macrophage"
genes <- c("TNFRSF1B","ITGA4","ITGB1","CD6","CXCL2" ,"CXCL3" , "CXCL8" , "GRN","SORT1","ITGB2" , "ICAM1",
           "ICAM2","NECTIN2" , "TIGIT","TNF" , "TNFRSF1A","TNF" , "VCAM1" )
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
# ythdf2 termcor----
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
selected_genes <- c("YTHDF2") 
selected_cell_type <- "Macrophage"
adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/GO:BP.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
selected_features <- c(  
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_MACROPHAGE_CYTOKINE_PRODUCTION",
  "GOBP_CYTOKINE_PRODUCTION",
  "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_CELL_CELL_ADHESION",
  "GOBP_REGULATION_OF_CELL_ADHESION",
  "GOBP_LEUKOCYTE_CELL_CELL_ADHESION",
  "GOBP_INFLAMMATORY_RESPONSE",
  "GOBP_LEUKOCYTE_CHEMOTAXIS_INVOLVED_IN_INFLAMMATORY_RESPONSE",
  "GOBP_MACROPHAGE_INFLAMMATORY_PROTEIN_1_ALPHA_PRODUCTION"
)
df <- filtered_adata %>%
  filter(feature_x %in% selected_genes,
         feature_y %in% selected_features) 
df$feature_y <- gsub("GOBP_", "", df$feature_y)
df$feature_y <- gsub("_", " ", df$feature_y)
df$feature_y <- tolower(df$feature_y)
p <- df %>%
  dplyr:: mutate(
    feature_x = fct_relevel(feature_x, selected_genes),
    feature_y = fct_relevel(feature_y, selected_features)
  ) %>%
  ggplot(
    aes(x = feature_x,
        y = feature_y)
  ) +
  geom_point(aes(size=-log10(p_value), color=estimate),
             shape = 19 
             # stroke = 3
  ) +
  scale_color_gradient2(low = "#0045d1", high = "red", mid = "yellow") +
  coord_fixed() +
  theme_bw()+
  theme(
    panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.text.x = element_text(face = "italic"),
    axis.line = element_line(color = "black",size = 0.2)
  ) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.ticks = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        plot.title = element_text(hjust = 0.5)) 
ggsave(
  "cor.pdf",
  plot = p,
  height = 6,
  width = 8,
)

