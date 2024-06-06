# compute-----
library(Seurat)
library(tidyverse)
source("custom_function.R")
seurat_obj_path <-
  "~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds"
selected_cell_type <- "Macrophage"
k <- 15
outdir <- "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/"
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "data")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
seurat_obj<- seurat_obj_path %>% readRDS()
seurat_obj <- seurat_obj
sub_seurat_obj <- seurat_obj %>%
  subset(cell_type == selected_cell_type)
sub_seurat_obj <- sub_seurat_obj %>%
  cat_construct_metacells(k = k, name = selected_cell_type)
metacell_obj <- hdWGCNA::GetMetacellObject(sub_seurat_obj)
table(metacell_obj$donor)
sub_seurat_obj %>% saveRDS(file.path(outdir, "seurat_obj.rds"))
head(sub_seurat_obj)

setwd("/home/pingxr/Atherosis_0723/code_human/figure/Macrophage/newdata_try/")
#### Load packages 
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
seurat_obj_path <-
  "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/data/seurat_obj.rds"
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2","IGF2BP3","YTHDF1","YTHDF2")
selected_cell_type <- "Macrophage"
k <- 15
outdir <- "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/"
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "correlation")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
seurat_obj <- seurat_obj_path %>% readRDS()
metacell_obj <- seurat_obj %>% hdWGCNA::GetMetacellObject()
metacell_obj %>% dim()
cell_types <- names(table(metacell_obj$cell_type))
selected_genes <-
  selected_genes[selected_genes %in% rownames(metacell_obj)]
all_genes <- rownames(metacell_obj)
genes_cor_res <- metacell_obj %>%
  cat_correlation(feature_x = selected_genes,
                  feature_y = all_genes)
write_csv(genes_cor_res,
          file = file.path(outdir,
                           "gene.cor.csv"))
category <- "H"   
pathway_cor_res <- metacell_obj %>%
  cat_correlation(feature_x = selected_genes,
                  feature_y = category,
                  outdir = outdir,  
  )  
write_csv(pathway_cor_res,
          file = file.path(outdir,
                           str_c(category, ".cor.csv")))
category <- "GO:BP"
pathway_cor_res1 <- metacell_obj %>%
  cat_correlation(feature_x = selected_genes,
                  feature_y = category,
                  outdir = outdir,
  )
write_csv(pathway_cor_res1,
          file = file.path(outdir,
                           str_c(category, ".cor.csv")))
category <- "CP:KEGG"
pathway_cor_res2 <- metacell_obj %>%
  cat_correlation(feature_x = selected_genes,
                  feature_y = category,
                  outdir = outdir,
  )
write_csv(pathway_cor_res2,
          file = file.path(outdir,
                           str_c(category, ".cor.csv")))

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")
selected_cell_type <- "Macrophage"
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2",
                    "IGF2BP3","YTHDF1","YTHDF2")
k <- 15
outdir <- "./Atherosis_0723/20230204_Atherosis/result/Macro"   
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "wgcna")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
seurat_obj <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/data/seurat_obj.rds")
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = selected_cell_type,
  group.by = "cell_type",
  assay = "RNA"
) 
seurat_obj <- TestSoftPowers(seurat_obj,
                             setDatExpr = FALSE)
plot_list <-
  PlotSoftPowers(seurat_obj)
p <- wrap_plots(plot_list, ncol = 2)
p
ggsave(
  file.path(outdir, "test_soft_powers.pdf"),
  plot = p,
  height = 6,
  width = 10
)
power_table <-
  GetPowerTable(seurat_obj)
head(power_table)
seurat_obj <-
  ConstructNetwork(seurat_obj,
                   soft_power = 10,  
                   setDatExpr = FALSE,
                   tom_outdir = outdir,
                   nThreads = 10)
pdf(file.path(outdir, "plot_dendrogram.pdf"))
PlotDendrogram(seurat_obj, main = str_c(selected_cell_type, " hdWGCNA Dendrogram"))
dev.off()
TOM <- GetTOM(seurat_obj)   
seurat_obj <-
  ModuleEigengenes(seurat_obj,
                   group.by.vars = "donor")
hMEs <-
  GetMEs(seurat_obj)
MEs <-
  GetMEs(seurat_obj, harmonized = FALSE)
seurat_obj <-
  ModuleConnectivity(seurat_obj,
                     group.by = "cell_type",
                     group_name = selected_cell_type)
p <-
  PlotKMEs(seurat_obj, ncol = 5)
p
ggsave(
  file.path(outdir, "plot_kmes.pdf"),
  plot = p,
  height = 10,
  width = 10
)
modules <-
  GetModules(seurat_obj)
modules %>%
  write_csv(file = file.path(outdir, "modules.csv"))
head(modules[, 1:6])
names(table(modules$module))
modules %>%
  dplyr::filter(gene_name %in% selected_genes) %>%
  group_by(module) %>%
  dplyr::count(gene_name) %>%
  write_csv(file = file.path(outdir, "modules_genes.csv"))
library(enrichR)
library(GeneOverlap)
dbs <- c("GO_Biological_Process_2021",
         "KEGG_2021_Human",
         "MSigDB_Hallmark_2020")
seurat_obj <-
  RunEnrichr(seurat_obj,
             dbs = dbs,
             
             max_genes = 100) 
saveRDS(seurat_obj, file = file.path(outdir, "wgcna.rds"))

enrichrtable <- hdWGCNA::GetEnrichrTable(seurat_obj)
enrichrtable %>%
  write_csv(file.path(outdir, "enrichr_table.csv"))

setwd("/home/pingxr/Atherosis_0723/20230204_Atherosis/code/Macro/")
library(Seurat)
library(tidyverse)
source("custom_function.R")
source("custom_function.R")
seurat_obj_path <-
  "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/data/seurat_obj.rds"
selected_genes <- c("METTL3")
selected_cell_type <- "Macrophage"
k <- 15
outdir <- "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/"
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "deg")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
seurat_obj <- seurat_obj_path %>% readRDS()
seurat_obj
DimPlot(seurat_obj,group.by = "cell_type")
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
deg_metacell_obj <- metacell_obj |> cat_deg(
  group_by = selected_genes,
  cell_types = selected_cell_type,
  min_pct = 0,
  mode = 'median',
  logfc_threshold = 0
)
table(metacell_obj@assays$RNA@counts['METTL3', ] > median(metacell_obj@assays$RNA@counts['METTL3', ]))
write_csv(deg_metacell_obj, file.path(outdir, "deg_metacell_obj.csv"))
enriched <- deg_metacell_obj |>
  filter(group == "METTL3") |>
  cat_enrich()
saveRDS(enriched, file.path(outdir, "METTL3_enriched.rds"))
gsea_res <- deg_metacell_obj |>
  filter(group == "METTL3") |>
  cat_gsea(category = "H"
  )
saveRDS(gsea_res, file.path(outdir, "h_gsea.rds"))
gsea_res1 <- deg_metacell_obj |>
  filter(group == "METTL3") |>
  cat_gsea(category = "C2"
  )
saveRDS(gsea_res1, file.path(outdir, "C2_gsea.rds"))
gsea_res2 <- deg_metacell_obj |>
  filter(group == "METTL3") |>
  cat_gsea(category = "C5"
  )
saveRDS(gsea_res2, file.path(outdir, "C5_gsea.rds"))

# pathway cor-----
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
# cor----
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
adata <-
  read_csv("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/correlation/gene.cor.csv")
filtered_adata <- adata %>%
  filter(p_value <= 0.05)
gene_sets <- msigdbr::msigdbr(species = "human",
                              category = "C5") %>%  
  dplyr::select(gs_name, gene_symbol)
selected_features <- gene_sets %>%
  filter(gs_name == "GOBP_MACROPHAGE_ACTIVATION") %>%
  pull(gene_symbol)

p <- filtered_adata %>%
  filter(feature_y %in% selected_features,
         feature_x == "METTL3") %>%
  mutate(change = if_else(
    p_value  <= 0.05 & abs(estimate) > 0,
    if_else(estimate > 0, "high", "low"),
    "Uncorrelation"
  )) %>%
  filter(change != "low") %>%
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
  size = 1,
  color = "#ff707a") +
  geom_point(aes(col = change), shape=18,size=3,colour="#0064b3",fill="#bde2ff") +
  theme_cat() +
  labs(title = "GOBP_MACROPHAGE_ACTIVATION",
       y = "Correlation of METTL3") +
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
  NoLegend() 
p   
ggsave(
  file.path(outdir, "METTL3-GOBP_MACROPHAGE_ACTIVATION.pdf"),
  plot = p,
  height = 6,
  width = 6
)
# volcano----
library('ggrastr')
library(ggplot2)
library(ggh4x)
library(cowplot)
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

seurat.obj <-readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/seurat/seurat_integration_anno.rds")
seurat.obj <- subset(seurat.obj,cell_type == 'Macrophage')
seurat.obj <- NormalizeData(seurat.obj)
expr <- as.matrix(seurat.obj@assays$RNA@data)
table <- AverageExpression(seurat.obj, group.by = "donor", assays = "RNA")[["RNA"]]
table <- as.data.frame(table)
a <- as.data.frame(rownames(table))
gene <- as.data.frame(c("METTL3", "RGS19","LPAR6","IGHM","IFNGR1","RB1","IGLC2",
                        "CMTM3","IFIT3","DEGS1","CX3CR1"))
colnames(gene) <- "name"
table1 <- table[c(gene$name),]
table1 <- t(table1)
corr <- round(cor(table1),3)
head(corr[,1:6])
p.mat <- round(cor_pmat(table1),3)
head(p.mat[,1:6])
p <- ggcorrplot(corr)
p
ggcorrplot(corr,method = "circle",type = "upper")
p <- ggcorrplot(corr, method = "circle",
                hc.order = TRUE,hc.method = "ward.D",
                outline.col = "white",ggtheme = theme_bw(),
                colors = c("#1442ff","white","#ff503f"),
                lab = TRUE,lab_size = 2
)
p
ggsave("cor2.pdf",
       p,height = 10,width = 10)

sce <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/Macro/Macrophage/k_15/wgcna/wgcna.rds")
pdf("tree.pdf",height = 4,width = 6)
PlotDendrogram(sce, main='INH hdWGCNA Dendrogram')
dev.off()
