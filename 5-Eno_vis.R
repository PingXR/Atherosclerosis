# cor----
library(Seurat)
library(tidyverse)
source("custom_function.R")
source("custom_plot_function.R")
setwd("~/20240103_Atherosis/v3/result/Fig_Eno/1-cor/2-m6Apathway")
outdir <- "~/20240103_Atherosis/v3/result/Fig_Eno/1-cor/2-m6Apathway/"
dir.create(outdir, recursive = T)
selected_genes <- c("METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
                    "YTHDC1","YTHDC2","YTHDF1","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
                    "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1") 
selected_cell_type <- "Endothelial"
adata <- read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/correlation/GO:BP.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
filtered_adata$feature_y <- gsub("GOBP_", "", filtered_adata$feature_y )
filtered_adata$feature_y  <- gsub("_", " ", filtered_adata$feature_y )
filtered_adata$feature_y  <- tolower(filtered_adata$feature_y )
selected_features <- c(     
  "cell migration involved in sprouting angiogenesis",
  "sprouting angiogenesis",
  "regulation of sprouting angiogenesis",
  "blood vessel endothelial cell proliferation involved in sprouting angiogenesis",
  "angiogenesis involved in wound healing",
  "blood vessel morphogenesis",
  "blood vessel endothelial cell migration",
  "endothelial cell migration",
  "positive regulation of endothelial cell proliferation"
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
  coord_fixed() +
  theme(
    panel.grid = element_line(size = 0.2, color = "lightgrey"),
    axis.text.x = element_text(face = "italic")
  )
p  
ggsave(
  file.path(outdir, "pathway_cor_dotplot.pdf"),
  plot = p,
  height = 4,
  width = 6,
)
#  cor2 ----
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
setwd("~/20240103_Atherosis/result/Fig_Eno/1-cor")
outdir <- "~/20240103_Atherosis/result/Fig_Eno/1-cor"

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
p <- expr %>%
  catscatter(
    x = ALKBH5,
    y = GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS,
    method = "pearson") +
  theme(aspect.ratio = 1)
ggsave("GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS1.pdf",height = 4,width = 4)
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
  ggplot(aes(x = ALKBH5, y = GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS)) +
  ggrastr::rasterise(  ggpointdensity::geom_pointdensity(size = 0.3),dpi = 600) +
  geom_smooth(          method = "lm",
                        formula = y ~ x,
                        color = "#005291",
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
  labs(y = "GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS" |>
         str_replace_all("GOBP_", "") |>
         str_replace_all("_", " ") |>
         str_to_sentence(),
       title = "R = 0.22, P = 2.4e-07")
p2
ggsave(
  file.path(outdir, "ALKBH5-GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS2.pdf"),
  plot = p2,
  height = 4,
  width = 4
)

# gene exsperssion -------
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
setwd("~/20240103_Atherosis/result/Fig_Eno/1-m6agene")
outdir <- "~/20240103_Atherosis/result/Fig_Eno/1-m6agene"

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
p <- VlnPlot(
  sub_seurat_obj,
  features = "NRP1", 
  group.by = "ALKBH5_group",
  assay = "RNA"
)
p
p <- p$data %>%
  mutate(ident = fct_relevel(ident, c("high", "low"))) %>%
  ggplot(aes(x = ident, y = NRP1, fill = ident)) +
  geom_violin(cex=0.5,color = "black")+
  geom_boxplot(width=0.1,cex=0.5,color = "black",fill = "white") +
  ggsignif::geom_signif(textsize = 2, comparisons = list(c("low", "high"))) +
  theme_cat() +
  labs(title = "ALKBH5",
       y = "NRP1") + theme(aspect.ratio = 2,
                           axis.title.x = element_blank()) +
  scale_fill_manual(values = c(low = "#2e86c1",
                               high = "#e32d32")) +
  NoLegend()
p
ggsave("ALKBH5-NRP1.pdf",height = 3,width = 3)

library(Seurat)
library(tidyverse)
source("custom_function.R")
source("custom_plot_function.R")
setwd("~/20240103_Atherosis/v2/result/Fig_Eno/1-cor/5-m6Agenescore")
outdir <- "~/20240103_Atherosis/v2/result/Fig_Eno/1-cor/5-m6Agenescore"
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
p <- VlnPlot(
  sub_seurat_obj,
  features = "KLF2", 
  group.by = "ALKBH5_group",
  assay = "RNA"
)
p
p <- p$data %>%
  mutate(ident = fct_relevel(ident, c("high", "low"))) %>%
  ggplot(aes(x = ident, y = KLF2, fill = ident)) +
  geom_violin(cex=0.5,color = "black")+
  geom_boxplot(width=0.1,cex=0.5,color = "black",fill = "white") +
  ggsignif::geom_signif(textsize = 2, comparisons = list(c("low", "high"))) +
  theme_cat() +
  labs(title = "ALKBH5",
       y = "KLF2") + theme(aspect.ratio = 2,
                           axis.title.x = element_blank()) +
  scale_fill_manual(values = c(low = "#78A08F",
                               high = "#DF935F")) +
  NoLegend()
p
ggsave("ALKBH5-KLF2.pdf",height = 3,width = 3)

# gene cor----
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
setwd("~/20240103_Atherosis/result/Fig_Eno/1-m6agene")
outdir <- "~/20240103_Atherosis/result/Fig_Eno/1-m6agene"

adata <-
  read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
gene_sets <- msigdbr::msigdbr(species = "human",
                              category = "C5") %>%  
  dplyr::select(gs_name, gene_symbol)
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS") %>%
  pull(gene_symbol)
p <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "ALKBH5") %>%
  mutate(change = if_else(
    p_value  <= 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation"
  )) %>%
  filter(change != "Uncorrelation") %>%
  top_n(40, wt = abs(estimate)) %>%
  mutate(feature_y = fct_reorder(feature_y, desc(estimate))) %>%
  ggplot(aes(x = feature_y, y = estimate)) +
  geom_segment(aes(
    x = feature_y,
    xend = feature_y,
    y = 0,   
    yend = estimate,
    color = change
  ),
  size = 1
  ) +
  geom_point(aes(col = change), shape=21,size=2,colour="#E32D32",fill="#EE7072") +
  theme_cat() +
  labs(title = "GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
       y = "Correlation of ALKBH5") +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(
      angle = 90,
      hjust = 0.5,
      face = "italic",
      vjust = 0.5
    ),
    axis.title.x = element_blank(),
    legend.position = c(0.55, 0.9),
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.margin = margin(b = -10)
  ) +
  geom_hline(yintercept = 0,lty=2,col="black",lwd=0.4) +
  scale_color_manual(values = c(
    low = "#0074a8",
    high = "#f19802"
  )) +
  NoLegend()
p   
ggsave(
  file.path(outdir, "ALKBH5-GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS.pdf"),
  plot = p,
  height = 4,
  width = 6
)

# volcano-----
library(Seurat)
library(tidyverse)
library(ggplot2)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
setwd("~/20240103_Atherosis/result/Fig_Eno/2-deg")
outdir <- "~/20240103_Atherosis/result/Fig_Eno/2-deg"

adata <- read_csv("~/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/deg/deg_metacell_obj.csv")
avg_log2FC=0.2
p_val = 0.05
deg_2 <- adata
k1 = (deg_2$p_val < p_val)&(deg_2$avg_log2FC < -avg_log2FC)
k2 = (deg_2$p_val < p_val)&(deg_2$avg_log2FC > avg_log2FC)
change = ifelse(k1,"down(16)",ifelse(k2,"up(20)","stable"))
table(change)
deg_2 <- deg_2[,colnames(deg_2)!="change"]
deg_2 <- mutate(deg_2,change=change)
table(deg_2$change)
dat  = deg_2
rownames(dat) <- dat$gene
for_label <- dat[c("ALKBH5","CLDN5","RNASE1","IL32","MYC","SOX17",
                   "POSTN","CD9","TAGLN","GAS6"),]
p <- ggplot(data = dat, 
            aes(x = avg_log2FC, 
                y = -log10(p_val))) +
  geom_point(alpha=0.5, size=2, 
             aes(color=change)) +
  ylab("-log10(p_val)")+
  scale_color_manual(values=c("#00AFBB", "#999999", "#FC4E07"))+
  geom_vline(xintercept=c(-avg_log2FC,avg_log2FC),lty=3,col="black",lwd=0.6) +
  geom_hline(yintercept = -log10(p_val),lty=3,col="black",lwd=0.6) +
  theme_bw() 
p
volcano_plot <- p +
  geom_point(size = 0.8, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = gene),
    data = for_label,
    size = 2,
    box.padding = 0.3,
    label.padding = 0.2,
    label.size = 0.1,
    label.r = 0.1,
    color="black"
  ) +
  theme_cat() +
  theme(
    aspect.ratio = 1,
    legend.position = "top",
    legend.title = element_blank(),
    legend.margin = margin(b = -8)
  )
volcano_plot
ggsave(filename = file.path(outdir, "volcano.pdf"),
       plot = volcano_plot,
       height = 4,
       width = 5)
# gsea----
library(Seurat)
library(tidyverse)
library(ggplot2)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
setwd("~/20240103_Atherosis/result/Fig_Eno/2-deg")
outdir <- "~/20240103_Atherosis/result/Fig_Eno/2-deg"

adata <- readRDS("~/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/deg/C5_gsea.rds")
a <-adata@result
a$Description <- gsub("GOBP_", "", a$Description)
a$Description <- gsub("_", " ", a$Description)
a$Description <- tolower(a$Description)
geneSetID = c('blood vessel endothelial cell proliferation involved in sprouting angiogenesis',
              'cell migration involved in sprouting angiogenesis',
              'sprouting angiogenesis',
              "positive regulation of blood vessel endothelial cell migration",
              'blood vessel endothelial cell migration',
              "positive regulation of vascular endothelial cell proliferation",
              "endothelial cell apoptotic process",
              "regulation of vasculogenesis"
)
dat <- subset(a,a$Description %in% geneSetID)
dat_plot <- data.frame(id = dat$Description,
                       nes = dat$NES,
                       GeneRatio = dat$setSize)
dat_plot$id <- factor(dat_plot$id,levels = dat_plot$id)
dat_plot$threshold <- ifelse(dat_plot$nes > 0,'Up','Down')
library(ggplot2)
theme = theme(axis.text.x = element_text(hjust = 0.5,colour = "black",size = 16), 
              axis.ticks.y = element_blank(), 
              axis.text.y = element_text(size = 18,face = "plain"), 
              axis.title.x = element_text(size = 18,face = "plain"), 
              axis.title.y = element_text(size = 18,face = "plain"), 
              axis.line = element_line(colour = "black",size = 1),
              plot.title = element_text(hjust = 0.5,size =  20),
              legend.title = element_blank(),
              legend.text = element_text(size = 20), 
              legend.position = c(0.9,0.15),
              legend.background = element_rect(fill = 'transparent'))
p <- ggplot(dat_plot, aes(x = reorder(id, order(nes, decreasing=F)), y = nes, fill = threshold)) + 
  geom_bar(stat = "identity") + 
  coord_flip() + 
  xlab('KEGG Pathway Enrichment') + 
  ylab(bquote("-"~Log[10]~"(P value)")) +
  geom_text(aes(label=GeneRatio),size =6) + 
  theme_bw() + theme
p
dat_plot <- dat_plot[order(dat_plot$nes,decreasing = T),]
p2 <- ggplot(data = dat_plot,aes(x = reorder(id, order(nes, decreasing=F)),y = nes,fill = threshold)) +
  geom_col()+
  coord_flip() +
  geom_text(aes(label = GeneRatio),size =5) + 
  scale_fill_manual(values = c('Up'= '#36638a','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-1,1),color = 'white',size = 0.5,lty='dashed') +
  xlab('GOBP Pathway Enrichment') + 
  ylab("NES") 
p2
up_pathway <- length(which(dat_plot$nes > 0))
down_pathway <- length(which(dat_plot$nes < 0))
high <- nrow(dat_plot)
p3 <- p2 + geom_text(data = dat_plot[1:up_pathway,],aes(x = id,y = 0.1,label = id),
                     hjust=1.2,color = 'black',size = 4) +
  geom_text(data = dat_plot[(up_pathway +1):high,],aes(x = id,y = -0.1,label = id),
            hjust = -0.2,color = 'black',size = 4) +
  scale_x_discrete(labels = NULL) +
  theme_bw() +theme(panel.grid=element_blank())
p3
ggsave("gsea2.pdf",p3,height = 6,width = 12)

# go-----
library(Seurat)
library(tidyverse)
library(ggplot2)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/tutorial/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
setwd("~/20240103_Atherosis/result/Fig_Eno/2-deg")
outdir <- "~/20240103_Atherosis/result/Fig_Eno/2-deg"

adata <- readRDS("~/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/deg/ALKBH5_enriched.rds")
filtered_adata <- adata %>%
  filter(P.value <= 0.05)
selected_features <- c(
  "positive regulation of endothelial cell development (GO:1901552)",
  "regulation of angiogenesis (GO:0045765)",
  'blood vessel morphogenesis (GO:0048514)',
  "Angiogenesis",
  "regulation of establishment of endothelial barrier (GO:1903140)",
  "cardiac endothelial cell differentiation (GO:0003348)",
  "positive regulation of wound healing (GO:0090303)"
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
    ylab = "ALKBH5 high cells enriched",
    hline = -log10(0.05),
    label_position = "inward"
  ) +
  scale_fill_manual(values = c(Downregulated = "#7ca878",
                               Upregulated = "#FC4E07")) +
  NoLegend()
p
ggsave(
  file.path(outdir, "enriched.pdf"),
  plot = p,
  height = 3,
  width = 4,
)

# hdwgcna-----
library(Seurat)
library(tidyverse)
library(igraph)
library(WGCNA)
library(ggplot2)
library(hdWGCNA)
library(AUCell)
setwd("~/20240103_Atherosis/result/Fig_Eno/3-wgcna/")
outdir <- c("~/20240103_Atherosis/result/Fig_Eno/3-wgcna/")
sce <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/wgcna/wgcna.rds")
pdf("tree.pdf",height = 4,width = 6)
PlotDendrogram(sce, main='INH hdWGCNA Dendrogram')
dev.off()

library(Seurat)
library(tidyverse)
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library(igraph)
library(WGCNA)
library(hdWGCNA)
setwd("~/20240103_Atherosis/result/Fig_Eno/3-wgcna")
outdir <- c("~/20240103_Atherosis/result/Fig_Eno/3-wgcna/")
seurat_obj <-
  readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/wgcna/wgcna.rds")
ModuleNetworkPlot(seurat_obj, outdir = file.path("./ModuleNetwork/"))
enrich_df <-
  GetEnrichrTable(seurat_obj)
enrich_df %>% filter(module == "turquoise", P.value <= 0.05) %>% View()
selected_features <- c(
  "Fluid shear stress and atherosclerosis",  
  "regulation of blood vessel endothelial cell migration (GO:0043535)",
  "regulation of angiogenesis (GO:0045765)",
  "positive regulation of endothelial cell migration (GO:0010595)",
  "regulation of cell migration involved in sprouting angiogenesis (GO:0090049)",
  "negative regulation of sprouting angiogenesis (GO:1903671)",
  "Lipid and atherosclerosis",
  "positive regulation of endothelial cell development (GO:1901552)",
  "positive regulation of blood vessel endothelial cell migration (GO:0043536)"
)
p <- enrich_df %>% filter(module == "turquoise",
                          P.value <= 0.05,
                          Term %in% selected_features) %>%
  mutate(
    change = "upregulated"
  ) %>%
  catbarplot(
    x = Term,
    y = -log10(P.value),
    sort = T,
    fill = change,
    group_by = change,
    flip = T,
    ylab = "Enriched",
    xlab = "-log10(P value)",
    hline = -log10(0.05),
    label_position = "inward"
  ) +
  scale_fill_manual(values = c(downregulated = "#00a6e1",
                               upregulated = "turquoise")) +
  NoLegend()
p
ggsave("turquoise_go.pdf",height = 4,width = 6)

library(Seurat)
library(tidyverse)
library(igraph)
library(WGCNA)
library(ggplot2)
library(hdWGCNA)
library(AUCell)
setwd("~/20240103_Atherosis/result/Fig_Eno/3-wgcna/corpathway")
outdir <- c("~/20240103_Atherosis/result/Fig_Eno/3-wgcna/corpathway/")
sce <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/wgcna/wgcna.rds")
C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)
a <- as.data.frame(unique(C5_gene_sets$gs_name))
selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% c(
    "GOBP_SPROUTING_ANGIOGENESIS",
    "GOBP_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
    "GOBP_REGULATION_OF_CELL_MIGRATION_INVOLVED_IN_SPROUTING_ANGIOGENESIS",
    "GOBP_ENDOTHELIAL_CELL_MIGRATION",
    "GOBP_ENDOTHELIAL_CELL_PROLIFERATION"
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
a <- sce@meta.data
b <- cbind(a,aucs)
sce@meta.data <- b
feature <- "ALKBH5"
sce <-
  AddModuleScore(sce, features = list(feature), name = "ALKBH5")
colnames(aucs)
cur_traits <- c( 
  "cell migration involved in sprouting angiogenesis",
  "endothelial cell migration",
  "endothelial cell proliferation",
  "regulation of cell migration involved in sprouting angiogenesis",
  "sprouting angiogenesis"
)
sce$"type" <- "Endothelial"
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
  high_color = '#fd7469',
  mid_color = '#f5e102',
  low_color = '#0291d7',
  plot_max = 0.2,
  combine=TRUE
)
ggsave("heatmap.pdf",p,height = 6,width = 8)

library(Seurat)
library(tidyverse)
library(patchwork)
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/computing/custom_function.R")
source("/home/pingxr/Atherosis_0723/code/m6a/new/20220802_scrna-m6A/code/visualization/custom_plot_function.R")
library(igraph)
library(WGCNA)
library(hdWGCNA)
setwd("~/20240103_Atherosis/result/Fig_Eno/3-wgcna")
outdir <- c("~/20240103_Atherosis/result/Fig_Eno/3-wgcna/")
seurat_obj <-
  readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Endo/Endothelial/k_20/wgcna/wgcna.rds")
pdf("hub30.pdf",height = 5,width = 6)
p <- HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 30, n_other=0,
  edge_prop = 0.75,
  mods = 'turquoise'
)
dev.off()
print(p)
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, 
  n_neighbors=15, 
  min_dist=0.1 
)
umap_df <- GetModuleUMAP(seurat_obj)
plot_list <- ModuleFeaturePlot(seurat_obj)
pdf("4-featureplot.pdf",height = 5, width =10)
p <- wrap_plots(plot_list, ncol=6)
print(p)
dev.off()

# tf----
library(tidyverse)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer) 
library(rstatix)
library(networkD3)
setwd("~/20240103_Atherosis/v2/result/Fig_Eno/4-tf")
outdir <- ("~/20240103_Atherosis/v2/result/Fig_Eno/4-tf")
adata <- read_csv( "/home/pingxr/Atherosis_0723/result_human/Endo/newdata_tryresult/Endothelial/pyscenic/2022_human-Atherosis_knn_20/Endothelial.tfs_targets.csv" )
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2","IGF2BP3","YTHDF1","YTHDF2")
head(adata)
df <- adata %>% filter(target %in% selected_genes)
df$m6a <- df$target
df$m6a[df$m6a == "METTL3"] <- "m6A_Writer"
df$m6a[df$m6a == c("METTL14")] <- "m6A_Writer"
df$m6a[df$m6a == c("WTAP") ] <- "m6A_Writer"
df$m6a[df$m6a == c("ALKBH5")] <- "m6A_Eraser"
df$m6a[df$m6a == c("FTO")] <- "m6A_Eraser"
df$m6a[df$m6a == c("IGF2BP2") ] <- "m6A_Reader"
df$m6a[df$m6a == c("IGF2BP3")] <- "m6A_Reader"
df$m6a[df$m6a == c("YTHDF1") ] <- "m6A_Reader"
df$m6a[df$m6a == c("YTHDF2")] <- "m6A_Reader"
df$value <- 1
df2<-aggregate(df[,4],by=list(df$tfs,df$target),FUN=sum)
df <- df[,-1]
names(df2)<-c("source","target","value")
names(df)<-c("source","target","value")
df2
df
data<-bind_rows(df2,df)
str(data)
Sankeylinks<-data
Sankeynodes<-data.frame(name=unique(c(Sankeylinks$source,Sankeylinks$target)))
Sankeynodes$index<-0:(nrow(Sankeynodes) - 1)
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="source",by.y="name")  
Sankeylinks<-merge(Sankeylinks,Sankeynodes,by.x="target",by.y="name")  
Sankeydata<-Sankeylinks[,c(4,5,3)] 
names(Sankeydata)<-c("Source","Target","Value")  
Sankeyname<-Sankeynodes[,1,drop=FALSE]  
p <- sankeyNetwork(Links=Sankeydata,Nodes=Sankeyname, Source ="Source",  
                   Target = "Target", Value = "Value", NodeID = "name",fontSize = 6, nodeWidth = 15,
                   nodePadding = 10)
saveNetwork(p,"sankey.html")



