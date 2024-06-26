# cellchat----
library(CellChat)
library(patchwork)
library(Seurat)
library(tidyverse)
source("custom_function.R")
setwd("~/20240103_Atherosis/v2/result/Fig_Eno/5-supply/cellchat2/")
outdir <- "~/20240103_Atherosis/v2/result/Fig_Eno/5-supply/cellchat2/"
seurat_obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
seurat_obj <- NormalizeData(seurat_obj)
selected_genes <- c("ALKBH5")
selected_cell_type <- "Endothelial"
k <- 20
expr <- seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[selected_genes, ])
seurat_obj[[str_c(selected_genes, "_group")]] <-
  if_else(expr[selected_genes,] > median(expr[selected_genes, ]),
          "high", "low")
print(table(seurat_obj[[str_c(selected_genes, "_group")]]))
group <- "high"
output.dir <-
  paste0(outdir, group, "/") # must have "/"
dir.create(output.dir, recursive = T)
seurat.obj <- seurat_obj[, seurat_obj[["ALKBH5_group"]] == group] #
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
group <- "low"
output.dir <-
  paste0(outdir, group, "/") # must have "/"
dir.create(output.dir, recursive = T)
seurat.obj <- seurat_obj[, seurat_obj[["ALKBH5_group"]] == group] #
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

# mutilcellchat----
cellchat.control  <- readRDS("~/20240103_Atherosis/v2/result/Fig_Eno/5-supply/cellchat2/high/cellchat-high.rds")
cellchat.case  <- readRDS("~/20240103_Atherosis/v2/result/Fig_Eno/5-supply/cellchat2/low/cellchat-low.rds")
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
gg1 <- netVisual_heatmap(cellchat,font.size = 5,font.size.title = 5,color.heatmap = c("#196774", "#EF6024"))
gg2 <- netVisual_heatmap(cellchat, measure = "weight",font.size = 5,font.size.title = 5,color.heatmap = c("#196774", "#EF6024"))
pdf("3-热图.pdf",width = 6,height = 3.5)
p <- gg1 + gg2
print(p)
dev.off()
gg1 <- rankNet(cellchat,
               mode = "comparison", stacked = T,
               do.stat = TRUE) +coord_flip() +scale_fill_manual(values = c('#88BBB0','#FFC885'))
gg2 <- rankNet(cellchat,
               mode = "comparison", stacked = F, 
               do.stat = TRUE) +coord_flip()+scale_fill_manual(values = c('#88BBB0','#FFC885'))
gg1+gg2
ggsave("5-堆叠柱形图.pdf", plot = gg1+gg2, width = 12,height = 10)
pdf("7-ALL.pdf",height = 15,width = 15)
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(2), targets.use = c(1,3:15),slot.name = "netP", title.name = paste0("Signaling pathways sending from EC - ", names(object.list)[i]), legend.pos.x = 10)
}
dev.off()
levels(object.list[[1]]@idents)
pathways.show <- c("ADGRE5")
weight.max <- getMaxWeight(object.list,slot.name = c('netP') ,attribute =pathways.show)
pdf("ADGRE5.pdf",height = 10,width = 10)
par(mfrow = c(1,2), xpd=TRUE)  
for (i in 1: length(object.list)) {
  netVisual_aggregate(object.list[[i]],signaling = pathways.show,layout = "chord",
                      edge.weight.max = weight.max[1],edge.width.max = 10,
                      signaling.name = paste(pathways.show,names(object.list)[i]),
                      pt.title = 4,title.space = 4,
  )
}
dev.off()
netVisual_bubble(cellchat,
                 sources.use = c(2), 
                 targets.use =  c(1,3:15),
                 signaling = c("CXCL","COLLAGEN","PECAM1","LAMININ","VEGF","ADGRE5"),
                 font.size = 6,font.size.title = 6,remove.isolate = T,
                 comparison = c(1:2), angle.x = 45) +  #comparison = c(1, 2)就是group里面有几组 2组 Before Late,
  scale_color_distiller(palette = "RdYlBu")
ggsave("8-气泡图5.pdf",height = 4,width = 6)
save(object.list, file = "~/20240103_Atherosis/v2/result/Fig_Eno/5-supply/cellchat2/cellchat_object.list.RData")
save(cellchat, file = "~/20240103_Atherosis/v2/result/Fig_Eno/5-supply/cellchat2/cellchat_merged_.RData")
# cor----
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2","YTHDF2")
seurat_obj <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/data/seurat_obj.rds")
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <-
  readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/correlation/GO:BP.scores.rds")
scores[1:2, 1:2]
if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]
a <- as.data.frame(colnames(expr))
p <- expr %>%
  catscatter(
    x = ALKBH5,
    y = GOBP_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY, 
    method = "pearson") +
  theme(aspect.ratio = 1)
ggsave(
  file.path(outdir, "ALKBH5-GOBP_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY.pdf"),
  plot = p,
  height = 5,
  width = 5
)
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2","YTHDF2")
seurat_obj <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/data/seurat_obj.rds")
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes,] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/correlation/GO:BP.scores.rds")
scores[1:2, 1:2]
metacell_obj[["C5"]] <- CreateAssayObject(scores)
if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]
p2 <- expr %>%
  ggplot(aes(x = ALKBH5, y = GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.3),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#0c4068",
                        fill = "#e5e5e5",
                        size = 0.5,
                        alpha = 0.5) +
  theme_cat() +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(face = "italic")) +
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_gradientn(colors = c("#074da3", "#068d0d","#f8a553","#fffc42", "#ff5046") )+
  labs(y = "GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       title = "R = 0.18, P = 1.2e-05")
p2
ggsave(
  file.path(outdir, "./2/ALKBH5-GOBP_VASCULAR_ENDOTHELIAL_GROWTH_FACTOR_RECEPTOR_SIGNALING_PATHWAY.pdf"),
  plot = p2,
  height = 4,
  width = 4
)
# gene-----
adata <-
  readRDS(
    "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/correlation/GO:BP.scores.rds"
  )
seurat_obj <-
  readRDS(
    "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/data/seurat_obj.rds"
  )
sub_seurat_obj <- seurat_obj %>%
  hdWGCNA::GetMetacellObject() %>%
  subset(cell_type == "Endothelial")
pacman::p_load(Seurat)
pacman::p_load(tidyverse)
adata[1:2, 1:2]
colnames(adata)
colnames(sub_seurat_obj)
sub_seurat_obj[["score"]] <-
  CreateAssayObject(counts = adata)
sub_seurat_obj
x <- "ALKBH5"
expr <- sub_seurat_obj %>%
  GetAssayData(assay = "RNA", slot = "data")
median(expr[x, ])
sub_seurat_obj[[str_c(x, "_group")]] <-
  if_else(expr[x,] > median(expr[x, ]),
          "high", "low")
print(table(sub_seurat_obj[[str_c(x, "_group")]]))
sub_seurat_obj
DefaultAssay(sub_seurat_obj) <- "score"
rownames(sub_seurat_obj)
library(ggplot2)
library(ggsignif) 
library(ggdist)
Custom.color <- c("#d33c54","#3fa2c6","#d6b55a")
Vec1 <- c("high", "low")
comb_list <- list()
for(i in 1:(length(Vec1)-1)) {
  for(j in (i+1):length(Vec1)) {
    comb <- combn(c(Vec1[i], Vec1[j]), 2)
    if(!any(comb[1,] == comb[2,])) {
      comb_list[length(comb_list)+1] <- list(comb)
    }
  }
}
load("~/20240103_Atherosis/v2/result/Fig_Eno/5-supply/cellchat2/cellchat_merged_.RData")
df.net <- subsetCommunication(cellchat, slot.name = 'net')
high <- df.net[[1]]
high <- subset(high,pathway_name == c("CXCL","COLLAGEN","PECAM1","LAMININ","VEGF","ADGRE5"))
ligand <- unique(high$ligand)
receptor <- unique(high$receptor)
genes <- unique(c(ligand,receptor))
genes <- c("PGF" ,"CXCL12","COL1A1","COL1A2","COL4A1","COL4A2","COL6A1","COL6A2",
           "COL6A3","COL9A3","LAMA2", "LAMA4", "LAMB1", "LAMB2", "LAMC1", "COL4A4",
           "LAMA5", "ADGRE5","PECAM1","FLT1","CXCR4", "CD44","SDC1","SDC4",  "CD55")
genes <- c("ITGA1","ITGB1","ITGA6","ITGA7","ITGA10","ITGAV","ITGB8","ITGB4")
for( i in genes){
  p <- VlnPlot(
    sub_seurat_obj,
    features = i, 
    group.by = "ALKBH5_group",
    assay = "RNA"
  ) 
  p
  data <- p$data %>%
    mutate(ident = fct_relevel(ident, c("high", "low")))
  colnames(data) <- c("gene","ident")
  p <- ggplot(data, aes(x = ident, y = gene,fill=ident)) +
    geom_jitter(mapping = aes(color=ident),width = .05, alpha = 0.5,size=0.9)+
    geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ 
    stat_halfeye(mapping = aes(fill=ident),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) +
    scale_color_manual(values = Custom.color)+  
    xlab("") +  
    ylab("") +   
    ggtitle(i)+  
    coord_flip()+
    theme(axis.ticks.x = element_line(size = 0,color = "black"),  
          panel.background = element_rect(fill = "white", color = "white"),  
          panel.grid.major.x = element_line(color = "gray", size = 0), 
          panel.grid.minor.y = element_blank(),
          panel.border = element_rect(color = "black", fill = NA,linewidth = 1),
          legend.position = "none",
          axis.title.x = element_text(size = 8), 
          axis.title.y = element_text(size = 8),
          axis.text.x = element_text(size = 8,hjust = 0.3), 
          axis.text.y = element_text(size = 8), 
          plot.title = element_text(hjust = 0.5)
    )+
    geom_signif(comparisons = comb_list,step_increase = .1,map_signif_level = TRUE,vjust = 0.5,hjust= 0,
    ) 
  p
  ggsave(paste0("./png/",i,".png"),p,height = 3.5,width = 3.5)
}
genes <- c("COL4A2","COL6A3","FLT1","LAMA5","ITGB1","ITGA1")
for( i in genes){
  p <- VlnPlot(
    sub_seurat_obj,
    features = i, 
    group.by = "ALKBH5_group",
    assay = "RNA"
  ) 
  p
  data <- p$data %>%
    mutate(ident = fct_relevel(ident, c("high", "low")))
  colnames(data) <- c("gene","ident")
  p <- ggplot(data, aes(x = ident, y = gene,fill=ident)) +
    geom_jitter(mapping = aes(color=ident),width = .05, alpha = 0.5,size=0.9)+
    geom_boxplot(position = position_nudge(x = 0.14),width=0.1,outlier.size = 0,outlier.alpha =0)+ 
    stat_halfeye(mapping = aes(fill=ident),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.6) + 
    scale_fill_manual(values = Custom.color)+   
    xlab("") +  
    ylab("") +  
    ggtitle(i)+ 
    coord_flip()+
    theme(axis.ticks.x = element_line(size = 0,color = "black"),  
          panel.background = element_rect(fill = "white", color = "white"),  
          panel.grid.major.x = element_line(color = "gray", size = 0), 
          panel.grid.minor.y = element_blank(), 
          panel.border = element_rect(color = "black", fill = NA,linewidth = 1), 
          legend.position = "none", 
          axis.title.x = element_text(size = 8),  
          axis.title.y = element_text(size = 8), 
          axis.text.x = element_text(size = 8,hjust = 0.3), 
          axis.text.y = element_text(size = 8), 
          plot.title = element_text(hjust = 0.5)
    )+
    geom_signif(comparisons = comb_list,step_increase = .1,map_signif_level = TRUE,vjust = 0.5,hjust= 0,
    ) 
  p
  ggsave(paste0("./pdf/",i,".pdf"),p,height = 3.5,width = 3.5)
}

