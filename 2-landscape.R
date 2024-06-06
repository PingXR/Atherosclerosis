# tsne--------
library(Seurat)
library(tidyverse)
library(tidydr)
library(ggrastr)
setwd("~/20240103_Atherosis/v2/result/Fig1")
path <-"~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2.rds"
seurat_obj <- readRDS(path)
DefaultAssay(seurat_obj) <- "RNA"

mytheme <- theme_void() + #空白主题，便于我们后期添加tSNE箭头
  theme(plot.margin = margin(5.5,15,5.5,5.5)) #画布空白页缘调整
a <- ggplot(tsne,aes(x= tSNE_1 , y = tSNE_2 ,color = cell_type,raster = TRUE) ) +  
  geom_point_rast(size = 0.4 , alpha =0.6 )+
  theme_dr(xlength = 0.2, #x轴长度
           ylength = 0.2, #y轴长度
           arrow = grid::arrow(length = unit(0.1, "inches"), #箭头大小/长度
                               ends = 'last', type = "closed")) + #箭头描述信息
  theme(panel.grid = element_blank())
label <- tsne %>%
  group_by(cell_type)%>%
  summarise(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2))
head(label)
b <- a +
  geom_text(data = label,
            aes(x = tSNE_1, y = tSNE_2, label = cell_type),
            color = 'black', size = 3.5)
b
mycol <- c("#B07AA1","#F89C74","#66C5CC","#75ACC3","#FE88B1","#80BA5A","#F6CF71","#EDC948","#B84D64",
           "#008695","#59A14F","#FF9DA7","#F28E2B","#DCB0F2","#EE7072","#E73F74","#f69896","#9B8E8C",
           "#E68310","#E15759")
c <- b +
  guides(color = guide_legend(override.aes = list(size = 3)))+ 
  scale_color_manual(values = mycol) +
  scale_fill_manual(values = mycol)
c
ggsave(
  file.path( "cell_type_tsne3.pdf"),
  plot = c,
  height =5,
  width = 6
)

DimPlot(seurat_obj,split.by = "donor",reduction = "tsne",ncol = 4,
        cols = c("#B07AA1","#F89C74","#66C5CC","#75ACC3","#FE88B1","#80BA5A","#F6CF71","#EDC948","#B84D64",
                 "#008695","#59A14F","#FF9DA7","#F28E2B","#DCB0F2","#EE7072","#E73F74","#f69896","#9B8E8C",
                 "#E68310","#E15759"),
        raster = F
)
ggsave("donortsne2.pdf",height = 5,width = 12)

# 基因表达量--------
library(Seurat)
library(tidyverse)
setwd("~/20240103_Atherosis/v2/result/Fig1")
path <-"~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2.rds"
seurat_obj <- readRDS(path)
table(seurat_obj$cell_type)
seurat_obj <- subset(seurat_obj,idents = c("Fibroblast 1","Endothelial","Macrophage","Fibromyocyte",
                                            "T cell","Smooth muscle cell","Pericyte 1","Pericyte 2","B cell",
                                            "Plasma cell 1","Fibroblast 2","Neuron","Plasma cell 2","NK cell",
                                            "Mast cell"))
seurat_obj
saveRDS(seurat_obj,"~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")

seurat.obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
expr <-
  AverageExpression(seurat.obj,
                    group.by = "cell_type", assays = "RNA")[["RNA"]]
features <- c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
              "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
              "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
seurat.obj$cellcode <- rownames(seurat.obj@meta.data)
strat_samp <- seurat.obj@meta.data %>% group_by(cell_type) %>% sample_frac(size = .03)
seurat_obj_sample <- seurat.obj[,strat_samp$cellcode]
seurat_obj_sample <- seurat_obj_sample[,order(seurat_obj_sample@meta.data$cell_type)]
annotation <- as.data.frame(seurat_obj_sample$cell_type)
colnames(annotation) <- "cell_type"
ann_colors=list()
annotation_colors = c("#B07AA1","#F89C74","#66C5CC","#75ACC3","#FE88B1","#80BA5A","#F6CF71","#EDC948","#B84D64",
                      "#FF9DA7","#DCB0F2","#EE7072","#E73F74","#f69896","#E68310")
names(annotation_colors) <- levels(factor(annotation$cell_type))
ann_colors[["cell_type"]] <- annotation_colors
anno_col <- data.frame(celltype=seurat_obj_sample$cell_type)
p <- pheatmap::pheatmap(
  expr[features, ],
  scale = "row",
  cellwidth = 10,
  cellheight = 10,
  fontsize = 7,
  annotation = annotation,
  annotation_colors = ann_colors,
  cluster_rows = F,
  cluster_cols = F,
  number_color = "white",
  color = colorRampPalette(c("#541798", "white", "#e2413f"))(50)
)
ggsave("m6Aexpression.pdf",p,height = 5,width = 4)

# 热图----
library(Seurat)
library(tidyverse)
setwd("~/20240103_Atherosis/v2/result/Fig1")
seurat_obj <- readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
DefaultAssay(seurat_obj) <- "RNA"
m6a_features <- c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
                  "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
                  "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
counts <- as.data.frame(seurat_obj[["RNA"]]@counts)
m6a_features  <- subset(m6a_features,subset = m6a_features %in% rownames(counts))
seurat_obj$cellcode <- rownames(seurat_obj@meta.data)
strat_samp <- seurat_obj@meta.data %>% group_by(cell_type) %>% sample_frac(size = .03)
seurat_obj_sample <- seurat_obj[,strat_samp$cellcode]
seurat_obj_sample <- seurat_obj_sample[,order(seurat_obj_sample@meta.data$cell_type)]
annotation <- as.data.frame(seurat_obj_sample$cell_type)
colnames(annotation) <- "cell_type"
ann_colors=list()
annotation_colors = c("#B07AA1","#F89C74","#66C5CC","#75ACC3","#FE88B1","#80BA5A","#F6CF71","#EDC948","#B84D64",
                      "#FF9DA7","#DCB0F2","#EE7072","#E73F74","#f69896","#E68310")
names(annotation_colors) <- levels(factor(annotation$cell_type))
ann_colors[["cell_type"]] <- annotation_colors
mat4=as.matrix(seurat_obj_sample[["RNA"]]@data[m6a_features,])
mat4=mat4/max(mat4)
anno_col <- data.frame(celltype=seurat_obj_sample$cell_type)
library(pheatmap)
bk=c(seq(0, 1, by = 0.01))
pheatmap(mat4,cluster_rows = F,cluster_cols = F,
         cellwidth = 0.6,
         cellheight = 5,
         fontsize = 4,
         show_colnames = F,
         show_rownames = T,
         annotation = annotation,
         annotation_colors = ann_colors,
         filename="m6a_heatmap_celltype2.pdf",
         color = c(
           colorRampPalette(colors = c("#FFFFFF", "#EE7072"))(length(bk))
         ),
         legend_breaks = seq(0, 1, 0.5),
         breaks = bk
)

# featureplot----
seurat.obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
library(viridis)
library(scCustomize)
features <- c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
              "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
              "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
p <- FeaturePlot(
  seurat.obj,
  features = features,
  order = T,
  ncol = 4,
  pt.size = 5,
  cols = brewer.pal(8, name = "OrRd"), #OrRd #Blues
  reduction = "tsne",
  raster = TRUE,
) +
  theme(aspect.ratio = 1)
ggsave("featureplot3.pdf",p,height = 18,width = 14)

# 基因表达量2热图----
library(Seurat)
library(tidyverse)
library(ggpubr)
setwd("~/20240103_Atherosis/v2/result/Fig1/")
seurat.obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
m6a_features <- c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
                  "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
                  "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
all_expr <-
  AverageExpression(seurat.obj,
                    assays = "RNA",
                    slot = "data",
                    group.by = "cell_type")[["RNA"]]
m6a_expr <-
  AverageExpression(
    seurat.obj,
    assays = "RNA",
    slot = "data",
    group.by = "cell_type",
    features = m6a_features
  )[["RNA"]]
bk <- c(seq(0.2, 0.59, by = 0.01), seq(0.6, 1, by = 0.01))
p <- cor(all_expr) %>%
  pheatmap::pheatmap(
    cellwidth = 6,
    cellheight = 6,
    border_color = "white",
    width = 0,
    treeheight_row = 0.5,
    treeheight_col = 0.5,
    fontsize = 4,
    color = c(
      colorRampPalette(colors = c('#945182', "yellow",'#e23b54'))(length(bk)) #FFCC66 ##1c4ebd '#541798', '#e23b54' # "#0f86a9", "white", "#ed8b10"
    ),
    legend_breaks = seq(0, 1, 0.2),
    breaks = bk,
    main = "All"
  )
ggsave("allgene1.pdf",p,height = 4,width = 4)
bk <- c(seq(0.8, 0.8999, by = 0.01), seq(0.9, 1, by = 0.01))
p <- cor(m6a_expr) %>%
  pheatmap::pheatmap(
    cellwidth = 6,
    cellheight = 6,
    border_color = "white",
    treeheight_row = 0.5,
    treeheight_col = 0.5,
    fontsize = 4,
    color = c(
      colorRampPalette(colors = c('#945182',"yellow" ,'#e23b54'))(length(bk)) #FFCC99  #1c4ebd
    ),
    legend_breaks = seq(0.8, 1, 0.05),
    breaks = bk,
    main = "m6A"
  )
ggsave("m6Agene2.pdf",p,height = 4,width = 4)

# upset----
sce <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
DefaultAssay(sce) <- "RNA"
Idents(sce) <- "cell_type"
markers <-
  FindAllMarkers(
    sce,
    assay = "RNA",
    only.pos = T,
    logfc.threshold = 0.5
  )
write_csv(markers, file = "~/20240103_Atherosis/v2/result/1-dealdata/figure/findallmaker_0.5_anno.csv")
markers <- read_csv("~/20240103_Atherosis/v2/result/1-dealdata/figure/findallmaker_0.5_anno.csv")
m6a_features <- c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
                  "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
                  "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
cell_markers_wider <- markers %>%
  select(cluster, gene) %>%
  bind_rows(tibble(cluster = "m6A", gene = m6a_features)) %>%
  mutate(n = 1) %>%
  pivot_wider(names_from = cluster, values_from = n) %>%
  replace(is.na(.), 0) %>%
  as.data.frame()
cell_markers_wider[1:2, 1:2]
pdf("upset1.pdf",height = 6,width = 7,onefile = F)
p <- UpSetR::upset(
  cell_markers_wider, #%>% 
  nsets = 20,
  matrix.color = "#FDC086", 
  main.bar.color = "#73BAD6",
  sets.bar.color = c("#EE7072","#66C5CC","#F89C74","#B07AA1","#FE88B1","#f69896","#F6CF71","#80BA5A","#EDC948","#75ACC3",
                     "#B84D64","#E68310","#DCB0F2","#E73F74","#FF9DA7","#998B95"),
  shade.color = "#F7A685",
  nintersects = 50,
  mb.ratio = c(0.6, 0.4),
  order.by = c("freq", "degree"),
  decreasing = c(TRUE, F),
  show.numbers = "no",
  point.size = 1,
  line.size = 0.5,
)
print(p)
dev.off()
# 花
library(plotrix)
setwd("~/20240103_Atherosis/v2/result/Fig1/supply")
flower_plot <- function(sample, otu_num, core_otu, start, a, b, r, ellipse_col, circle_col) {
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type='n')
  n   <- length(sample)
  deg <- 360 / n
  res <- lapply(1:n, function(t){
    draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180), 
                 y = 5 + sin((start + deg * (t - 1)) * pi / 180), 
                 col = ellipse_col[t],
                 border = ellipse_col[t],
                 a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         otu_num[t])
    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = 1
      )
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = 1
      )
    }
  })
  draw.circle(x = 5, y = 5, r = r, col = circle_col, border = NA)
  text(x = 5, y = 5, paste('Core:', core_otu))
}
sce <- readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
v1 <- c("Fibroblast 1","Endothelial","Macrophage","Fibromyocyte","T cell","Smooth muscle cell","Pericyte 1",
        "Pericyte 2", "B cell","Plasma cell 1","Fibroblast 2", "Neuron","Plasma cell 2",
        "NK cell","Mast cell","m6A")
v2 <- c(174,288,307,66,97,87,111,71,77,17,74,309,20,96,120,23)
v3 <- 0
pdf("flower.pdf",height = 6,width = 6)
p <- flower_plot(sample = v1, 
                 otu_num = v2, 
                 core_otu = 0, 
                 start = 90, a = 0.5, b = 2, r = 1, 
                 ellipse_col = c("#B07AA1","#F89C74","#66C5CC","#75ACC3","#FE88B1","#80BA5A","#F6CF71","#EDC948","#B84D64",
                                 "#FF9DA7","#DCB0F2","#EE7072","#E73F74","#f69896",
                                 "#E68310","#998B95"), 
                 circle_col = 'white')
print(p)
dev.off()

# 气泡图----
seurat <- readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2.rds")
features <- c(
  "COL1A1","LUM","OMD", 
  "PECAM1","CLDN5",
  "CD163","CD14", 
  "GAS6","TNFRSF11B",
  "CD3E","TRAC",
  "ACTA2","TAGLN",
  "TRPC6","PDGFRB","RERGL","NET1", 
  "MS4A1","CD79B",            
  "APOE",           
  "S100A8",          
  "SLAMF7",  "TXNDC5","EAF2", "IGHG1",
  "S100B","PLP1",
  "KLRD1","GNLY", 
  "TPSAB1","CPA3",
  "EDN1","GZMB","TFF3"  
)
p4 <- DotPlot(seurat, features = features,
              col.min = 0) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = -45,hjust = 1,vjust = 1,size = 6),
        axis.text.y = element_text(size = 6)) +
  scale_color_gradientn(colours = c('#864b76',"yellow" ,'#e23b54')) 
pdf(file = "p4.pdf",width =10,height = 5)
print(p4)
dev.off()

# 比例----
seurat.obj <-
  readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
DefaultAssay(seurat.obj) <- "RNA"
seurat.obj <- NormalizeData(seurat.obj)
features <- c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
              "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
              "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
features[!(features %in% colnames(seurat.obj@assays$RNA@data))]
prop <- data.frame()
for (feature in features) {
  if (feature %in% rownames(seurat.obj@assays$RNA@data)) {
    pn <-
      as.data.frame(t(ifelse(
        as.matrix(seurat.obj[feature, ]@assays$RNA@data) > 0,
        "postive",
        "negative"
      )))
    seurat.obj <-
      AddMetaData(seurat.obj,
                  metadata = pn,
                  col.name = paste0("pn_", feature))
    df <-
      as.data.frame(prop.table(
        table(seurat.obj$cell_type,
              seurat.obj@meta.data[, paste0("pn_", feature)]),
        margin = 1
      ))
    df$feature <- feature
    prop <- rbind(prop, df)
  }
}
head(prop)
colnames(prop) <- c("cell_type", "group", "freq", "feature")
library(showtext)
pdf("percent.pdf",height = 3.5,width = 7)
p1 <-
  ggplot(prop[prop$group == "postive", ],
         aes(x = freq * 100, y = feature, fill = cell_type)
  ) +
  geom_bar(stat = "identity",color = "black",linewidth =0.1) +
  facet_grid(. ~ cell_type) +
  theme_bw() +
  theme(
    axis.text.y = element_text(
      face = "italic",
      colour = "black",
      size = 3,
    ),
    axis.title.y = element_blank(),
    axis.text.x = element_text(
      colour = "black",
      size = 3,
      # family = "Arial"
    ),
    strip.background = element_blank(),
    strip.text = element_text(size = 5, ),
    axis.title.x = element_text(size = 5, ),
    panel.grid = element_blank()
  ) +
  labs(x = "Cell proportion (%)") +
  guides(fill = "none") +
  scale_fill_manual(values = rep(c("#C7E9C0","#A1D99B","#00CC33","#74C476","#41AB5D","#238B45","#33CC33","#66CC00","#33CC00","#66CC66","#66CC99","#339966","#339900","#669900","#009966","#339933","#339900"),2))
print(p1)
dev.off()



