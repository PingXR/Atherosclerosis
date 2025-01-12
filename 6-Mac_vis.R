# cor-----
library(Seurat)
library(tidyverse)  
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library(tidyverse)
outdir <- "~/20240103_Atherosis/result/Fig_Mac/1-cor/"
dir.create(outdir, recursive = T)
setwd("~/20240103_Atherosis/result/Fig_Mac/1-cor")
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2","IGF2BP3","YTHDF1","YTHDF2") 
selected_cell_type <- "Macrophage"

selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2","IGF2BP3","YTHDF1","YTHDF2") 
selected_cell_type <- "Macrophage"
adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/GO:BP.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
selected_features <- c(    
  "GOBP_REGULATION_OF_MACROPHAGE_DIFFERENTIATION",
  "GOBP_MACROPHAGE_DIFFERENTIATION",
  "GOBP_REGULATION_OF_MACROPHAGE_PROLIFERATION",
  "GOBP_MACROPHAGE_CHEMOTAXIS",
  "GOBP_MACROPHAGE_ACTIVATION",
  "GOBP_REGULATION_OF_MACROPHAGE_DERIVED_FOAM_CELL_DIFFERENTIATION",
  "GOBP_MACROPHAGE_MIGRATION",
  "GOBP_LEUKOCYTE_DEGRANULATION",
  "GOBP_IMMUNE_RESPONSE",
  "GOBP_REGULATION_OF_INFLAMMATORY_RESPONSE",
  "GOBP_LEUKOCYTE_TETHERING_OR_ROLLING"
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
p  
ggsave(
  file.path(outdir, "pathway_cor_dotplot-gobp2.pdf"),
  plot = p,
  height = 6,
  width = 8,
)

# pathway cor----
library(Seurat)
library(tidyverse) 
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library(tidyverse)
outdir <- "~/20240103_Atherosis/result/Fig_Mac/1-cor/"
setwd("~/20240103_Atherosis/result/Fig_Mac/1-cor")
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2","IGF2BP3","YTHDF1","YTHDF2") 
selected_cell_type <- "Macrophage"

seurat_obj <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/data/seurat_obj.rds")
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
expr <-
  GetAssayData(metacell_obj, slot = "data", assay = "RNA")[selected_genes, ] %>%
  as.data.frame()
expr[1:2, 1:2]
scores <-
  readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/GO:BP.scores.rds")
scores[1:2, 1:2]
if (identical(colnames(scores), colnames(expr))) {
  expr <- rbind(scores, expr)
  expr <- expr %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
}
expr[1:2, 1:2]
p <- expr %>%
  catscatter(
    x = METTL3,  
    y = GOBP_REGULATION_OF_MACROPHAGE_DIFFERENTIATION,
    method = "pearson") +
  theme(aspect.ratio = 1)
p
ggsave(
  file.path(outdir, "METTL3-GOBP_REGULATION_OF_MACROPHAGE_DIFFERENTIATION.pdf"),
  plot = p,
  height = 5,
  width = 5
)
library(rcartocolor)
p2 <- expr %>%
  ggplot(aes(x = METTL3, y = GOBP_REGULATION_OF_MACROPHAGE_DIFFERENTIATION)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#4552A0",
                        fill = "#b9d6e9",
                        size = 0.5,
                        alpha = 0.3) +
  theme_cat() +
  theme(aspect.ratio = 1, legend.margin = margin(l=-8),
        axis.title.x = element_text(face = "italic")) +
  guides(color = guide_colorbar(
    frame.colour = "black",
    frame.linewidth = 0.5,
    ticks.colour = "black",
    title = "Density"))+
  scale_color_continuous(low = "#0035e6", high = "#ff5046")+
  labs(y = "GOBP_REGULATION_OF_MACROPHAGE_DIFFERENTIATION" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       title = "R = 0.32, P = 7.8e-19")
p2
ggsave(
  file.path(outdir, "METTL3-GOBP_REGULATION_OF_MACROPHAGE_DIFFERENTIATION2.pdf"),
  plot = p2,
  height = 5,
  width = 4,
)

# genecor----
library(Seurat)
library(tidyverse) 
library(ggridges)
outdir <- "~/20240103_Atherosis/v3/result/Fig_Mac/1-cor/4-genescore"
setwd("~/20240103_Atherosis/v3/result/Fig_Mac/1-cor/4-genescore")

adata <-
  readRDS(
    "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/GO:BP.scores.rds"
  )
seurat_obj <-
  readRDS(
    "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/data/seurat_obj.rds"
  )
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
rownames(sub_seurat_obj)
genes <- c("TRIB1","RB1","IFNGR1","CX3CR1")
df <- data_frame()
for(i in genes){
  p <- VlnPlot(
    sub_seurat_obj,
    features = i,
    group.by = "METTL3_group",
    assay = "RNA"
  )
  data <- p$data
  data$gene <- paste0(i)
  colnames(data) <- c("value","group","gene")
  df <- rbind(df,data)
}
table(df$gene)
table(df$group)
cols=c("#EE757A","#78ADD5","#A78CC1", "#DBCA85","#FFAE5D","#CE9D7D","#84C780","#DED2E6","#C8E2EF")
p <- ggplot(df, aes(x = value, y = gene,fill=group,color=group)) + 
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
  theme(axis.title.x = element_text(size = 5,color = "black"),  
        axis.title.y = element_text(size = 5,color = "black"),  
        axis.text.x = element_text(size = 5,color = "black"),   
        axis.text.y = element_text(size = 5,color = "black"))+
  theme(plot.title = element_text(hjust = 0.2, vjust = 2,size = 5)  
  )
ggsave("./supply/mettl3.pdf",p,height = 4.5,width = 4.5)

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
rownames(sub_seurat_obj)
genes <- c("TRIB1","RB1","IFNGR1","CX3CR1")
df <- data_frame()
for(i in genes){
  p <- VlnPlot(
    sub_seurat_obj,
    features = i,
    group.by = "YTHDF2_group",
    assay = "RNA"
  )
  data <- p$data
  data$gene <- paste0(i)
  colnames(data) <- c("value","group","gene")
  df <- rbind(df,data)
}
table(df$gene)
table(df$group)
cols=c("#FFAE5D","#A78CC1", "#DBCA85","#78ADD5","#CE9D7D","#84C780","#DED2E6","#C8E2EF")
p <- ggplot(df, aes(x = value, y = gene,fill=group,color=group)) + 
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
  theme(axis.title.x = element_text(size = 5,color = "black"),  
        axis.title.y = element_text(size = 5,color = "black"),  
        axis.text.x = element_text(size = 5,color = "black"),   
        axis.text.y = element_text(size = 5,color = "black"))+
  theme(plot.title = element_text(hjust = 0.2, vjust = 2,size = 5)  
  )
ggsave("./supply/ythdf2.pdf",p,height = 4.5,width = 4.5)

# volcano----
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library('ggrastr')
library(ggplot2)
library(ggh4x)
library(cowplot)
outdir <- "~/20240103_Atherosis/result/Fig_Mac/2-deg"
setwd("~/20240103_Atherosis/result/Fig_Mac/2-deg")

adata <- read.csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/deg-METTL3/deg_metacell_obj.csv")
avg_log2FC=0.2
p_val = 0.05
deg_2 <- adata
k1 = (deg_2$p_val < p_val)&(deg_2$avg_log2FC < -avg_log2FC)
k2 = (deg_2$p_val < p_val)&(deg_2$avg_log2FC > avg_log2FC)
change = ifelse(k1,"Downregulated(124)",ifelse(k2,"Upregulated(142)","stable"))
table(change)
deg_2 <- deg_2[,colnames(deg_2)!="change"]
deg_2 <- mutate(deg_2,change=change)
table(deg_2$change)
dat  = deg_2
rownames(dat) <- dat$gene
for_label <- dat[c("METTL3","TMSB4X","PPIC","PABPC1","MYL6","MIF","EIF5A",
                   "SERF2","EIF1","CD63","GAS5","S100A10"),]
p <- ggplot(data = dat,
            aes(x = avg_log2FC,
                y = -log10(p_val))) +
  geom_point_rast(  alpha=0.5, 
                    aes(color=change),
                    size=2,
                    raster.dpi = getOption("ggrastr.default.dpi", 300))+
  ylab("-log10(p_val)")+
  scale_color_manual(values=c("#0064B3", "#999999", "#f84200"))+
  geom_vline(xintercept=c(-avg_log2FC,avg_log2FC),lty=3,col="black",lwd=0.6) +
  geom_hline(yintercept = -log10(p_val),lty=3,col="black",lwd=0.6) +
  theme_bw()
p
volcano_plot <- p +
  geom_point(size = 0.1, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = gene),
    data = for_label,
    size = 2,
    box.padding = 0.3,
    label.padding = 0.2,
    label.size = 0.1,
    label.r = 0.1,
    color="black",
  ) +
  theme_cat() +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    legend.title = element_blank(),
    legend.margin = margin(b = -8)
  )
volcano_plot
ggsave(filename = file.path(outdir, "METTL3volcano2.pdf"),
       plot = volcano_plot,
       height = 5,
       width = 5)
adata <- read.csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/deg-YTHDF2//deg_metacell_obj.csv")
avg_log2FC=0.2
p_val = 0.05
deg_2 <- adata
k1 = (deg_2$p_val < p_val)&(deg_2$avg_log2FC < -avg_log2FC)
k2 = (deg_2$p_val < p_val)&(deg_2$avg_log2FC > avg_log2FC)
change = ifelse(k1,"Downregulated(270)",ifelse(k2,"Upregulated(299)","stable"))
table(change)
deg_2 <- deg_2[,colnames(deg_2)!="change"]
deg_2 <- mutate(deg_2,change=change)
table(deg_2$change)
dat  = deg_2
rownames(dat) <- dat$gene
for_label <- dat[c("YTHDF2","HCLS1","ID2","TRIB1","RB1","PTPN2",
                   "CD63","GAS5","S100A10"),]
p <- ggplot(data = dat, 
            aes(x = avg_log2FC, 
                y = -log10(p_val)))+
  geom_point_rast(  alpha=0.5, 
                    aes(color=change),
                    size=2,
                    raster.dpi = getOption("ggrastr.default.dpi", 300))+
  ylab("-log10(p_val)")+
  scale_color_manual(values=c("#0064B3", "#999999", "#f84200"))+
  geom_vline(xintercept=c(-avg_log2FC,avg_log2FC),lty=3,col="black",lwd=0.6) +
  geom_hline(yintercept = -log10(p_val),lty=3,col="black",lwd=0.6) +
  theme_bw()
p
volcano_plot <- p +
  geom_point(size = 0.1, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = gene),
    data = for_label,
    size = 2,
    box.padding = 0.3,
    label.padding = 0.2,
    label.size = 0.1,
    label.r = 0.1,
    color="black",
  ) +
  theme_cat() +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    legend.title = element_blank(),
    legend.margin = margin(b = -8)
  )
volcano_plot
ggsave(filename = file.path(outdir, "YTHDF2-volcano2.pdf"),
       plot = volcano_plot,
       height = 5,
       width = 5)
# go-----
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
outdir <- "~/20240103_Atherosis/result/Fig_Mac/2-deg"
setwd("~/20240103_Atherosis/result/Fig_Mac/2-deg")

adata_METTL3 <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/deg-METTL3/METTL3_enriched.rds")
filtered_adata <- adata_METTL3 %>%
  filter(P.value <= 0.05)
selected_features <- c(
  "macrophage activation (GO:0042116)",
  "regulation of myeloid cell differentiation (GO:0045637)",
  "positive regulation of macrophage differentiation (GO:0045651)",
  "regulation of cell migration (GO:0030334)",
  "regulation of cell cycle (GO:0051726)",
  "cellular response to type I interferon (GO:0071357)"
)
p <- filtered_adata %>%
  filter(Term %in% selected_features) %>%
  catbarplot(
    x = Term,
    y = -log10(P.value),
    sort = T,
    fill = change,
    group_by = change,
    flip = T,
    ylab = "METTL3 high cells enriched",
    hline = -log10(0.05),
    label_position = "inward"
  ) +
  scale_fill_manual(values = c(Downregulated = "#b9ef97",
                               Upregulated = "#ffbc95")) +
  NoLegend()
p
ggsave("METTL3-enriched.pdf",p,height = 4,width = 5)

adata_YTHDF2 <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/deg-YTHDF2/YTHDF2_enriched.rds")
filtered_adata <- adata_YTHDF2 %>%
  filter(P.value <= 0.05)
selected_features <- c(
  "positive regulation of monocyte differentiation (GO:0045657)",
  "regulation of cell differentiation (GO:0045595)",
  "macrophage activation (GO:0042116)",
  "macrophage differentiation (GO:0030225)",
  "macrophage migration (GO:1905517)",
  "regulation of macrophage proliferation (GO:0120041)",
  "cellular response to macrophage colony-stimulating factor stimulus (GO:0036006)",
  "response to macrophage colony-stimulating factor (GO:0036005)",
  "regulation of macrophage chemotaxis (GO:0010758)",
  "TNF-alpha Signaling via NF-kB",
  "TNF signaling pathway"
)
p <- filtered_adata %>%
  filter(Term %in% selected_features) %>%
  catbarplot(
    x = Term,
    y = -log10(P.value),
    sort = T,
    fill = change,
    group_by = change,
    flip = T,
    ylab = "YTHDF2 high cells enriched",
    hline = -log10(0.05),
    label_position = "inward"
  ) +
  scale_fill_manual(values = c(Downregulated = "#009adb",
                               Upregulated = "#ffb1c0")) +
  NoLegend()
p
ggsave("YTHDF2-enriched.pdf",p,height = 4,width = 6)

# gsea----
library(Seurat)
library(tidyverse)
library(GseaVis)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
outdir <- "~/20240103_Atherosis/result/Fig_Mac/2-deg"
setwd("~/20240103_Atherosis/result/Fig_Mac/2-deg")

adata <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/deg-METTL3/C5_gsea.rds")  #不是一个表格，是一个gsea对象 class(adata)
View(adata@result)
geneSetID = c(
  "GOBP_MACROPHAGE_ACTIVATION",
  "GOBP_MACROPHAGE_PROLIFERATION",
  "GOBP_MACROPHAGE_DIFFERENTIATION")
gene <- c("TNF","CSF1R","TRIB1","CX3CR1","JUN","CCL3","THBS1")
p <- gseaNb(object = adata,
            geneSetID = geneSetID,
            addGene = gene,
            subPlot = 2,
            termWidth = 35,
            addPval = T,
            legend.position = c(0.8,0.8),
            pvalX = 0.08,
            pvalY = 0.2,
            curveCol = c("#b9770e", "#EB4747", "#35a132"),
            htCol = c("#2e86c1","palevioletred1"),
            segCol = "gray90"
)
p
ggsave(
  file.path(outdir, "gsea_METTL3.pdf"),
  plot = p,
  height = 8,
  width = 10)

adata <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/deg-YTHDF2/C5_gsea.rds")  #不是一个表格，是一个gsea对象 class(adata)
View(adata@result)
geneSetID = c(
  "GOBP_MACROPHAGE_ACTIVATION",
  "GOBP_MACROPHAGE_PROLIFERATION",
  "GOBP_MACROPHAGE_DIFFERENTIATION")
gene <- c("TNF","CSF1R","TRIB1","CX3CR1","JUN","CCL3","THBS1")
p <- gseaNb(object = adata,
            geneSetID = geneSetID,
            addGene = gene,
            subPlot = 2,
            termWidth = 35,
            addPval = T,
            legend.position = c(0.8,0.8),
            pvalX = 0.08,
            pvalY = 0.2,
            curveCol = c("#b9770e", "#EB4747", "#35a132"),
            htCol = c("#2e86c1","palevioletred1"),
            segCol = "gray90"
)
p
ggsave(
  file.path(outdir, "gsea_YTHDF2.pdf"),
  plot = p,
  height = 8,
  width = 10)

# hdwgcna-----
library(Seurat)
library(tidyverse)
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library(igraph)
library(WGCNA)
library(hdWGCNA)
setwd("~/20240103_Atherosis/result/Fig_Mac/3-wgcna")
outdir <- c("~/20240103_Atherosis/result/Fig_Mac/3-wgcna")

seurat_obj <-
  readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/wgcna/wgcna.rds")
ModuleNetworkPlot(seurat_obj, outdir = file.path("./ModuleNetwork/"))
enrich_df <-
  GetEnrichrTable(seurat_obj)
enrich_df %>% filter(module == "green", P.value <= 0.05) %>% View()
selected_features <- c(
  "regulation of inflammatory response (GO:0050727)",
  "cytokine-mediated signaling pathway (GO:0019221)",
  "myeloid cell differentiation (GO:0030099)",
  "positive regulation of cell activation (GO:0050867)",
  "positive regulation of myeloid leukocyte mediated immunity (GO:0002888)",
  "positive regulation of interferon-gamma production (GO:0032729)",
  "cellular response to cytokine stimulus (GO:0071345)"
)
p <- enrich_df %>% filter(module == "green",
                          P.value <= 0.05,
                          Term %in% selected_features) %>%
  catbarplot(
    x = Term,
    y = -log10(P.value),
    sort = T,
    fill = "green",
    flip = T,
    ylab = "Enriched",
    xlab = "-log10(P value)",
    hline = -log10(0.05),
    label_position = "inward"
  ) +
  NoLegend() 
p
ggsave(
  file.path(outdir, "METTL3-green_enriched.pdf"),
  plot = p,
  height = 4,
  width = 6
)

library(Seurat)
library(tidyverse)
library(igraph)
library(WGCNA)
library(ggplot2)
library(hdWGCNA)
library(AUCell)
setwd("~/20240103_Atherosis/result/Fig_Mac/3-wgcna/score")
outdir <- c("~/20240103_Atherosis/result/Fig_Mac/3-wgcna/score/")

sce <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/wgcna/wgcna.rds")
C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)
a <- as.data.frame(unique(C5_gene_sets$gs_name))
selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% c(
    "GOBP_MACROPHAGE_ACTIVATION",
    "GOBP_MACROPHAGE_DIFFERENTIATION",
    "GOBP_MACROPHAGE_PROLIFERATION",
    "GOBP_INFLAMMATORY_RESPONSE",
    "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE"
  )) 
selected_gene_sets
cells_rankings <- AUCell_buildRankings(as.matrix(sce@assays$RNA@data))
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
colnames(aucs) <- gsub("GOBP_", "", colnames(aucs))
colnames(aucs) <- gsub("_", " ", colnames(aucs))
colnames(aucs) <- tolower(colnames(aucs))
colnames(aucs)
a <- sce@meta.data
b <- cbind(a,aucs)
sce@meta.data <- b
feature <- "METTL3"
sce <-
  AddModuleScore(sce, features = list(feature), name = "METTL3")
cur_traits <- c( 
  "cytokine production involved in inflammatory response" ,                     
  "inflammatory response",
  "macrophage activation",                               
  "macrophage differentiation" ,        
  "macrophage proliferation"   ,
  "METTL31"
)
sce$"type" <- "Macrophage"
seurat_obj <- ModuleTraitCorrelation(
  sce,
  traits = cur_traits,
  group.by='type'
)
p <- PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = '#ff8d81',
  mid_color = '#ffec8d',
  low_color = '#7db2eb',
  plot_max = 0.2,
  combine=TRUE
)
ggsave("heatmap2.pdf",p,height = 6,width = 8)





