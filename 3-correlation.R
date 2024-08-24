# 热图-----
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(psych)
library(reshape2)
setwd("~/20240103_Atherosis/v2/result/Fig2/heatmap/")
sce <- readRDS("~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2sub.rds")
rt <- as.matrix(sce@assays$RNA@data)
m6a_features <- c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
                  "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
                  "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%
  dplyr::select(gs_name, gene_symbol)
selected_gene_sets <- C5_gene_sets %>%
  filter(gs_name %in% c(
    "GOBP_INFLAMMATORY_RESPONSE",
    "GOBP_ACUTE_INFLAMMATORY_RESPONSE",
    "GOBP_CHRONIC_INFLAMMATORY_RESPONSE",
    "GOBP_LIPID_BIOSYNTHETIC_PROCESS",
    "GOBP_LIPID_CATABOLIC_PROCESS",
    "GOBP_LIPID_OXIDATION",
    "GOBP_OXIDATIVE_PHOSPHORYLATION",
    "GOBP_CELL_PROLIFERATION_INVOLVED_IN_HEART_MORPHOGENESIS",
    "GOBP_CELL_MIGRATION_INVOLVED_IN_HEART_DEVELOPMENT",
    "GOBP_SPROUTING_ANGIOGENESIS",
    "GOBP_APOPTOTIC_SIGNALING_PATHWAY",
    "GOBP_RESPONSE_TO_FLUID_SHEAR_STRESS",
    "GOBP_CELL_CELL_ADHESION",
    "GOBP_CELL_ADHESION_INVOLVED_IN_HEART_MORPHOGENESIS",
    "GOBP_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY_INVOLVED_IN_HEART_DEVELOPMENT",
    "GOBP_ANGIOGENESIS_INVOLVED_IN_CORONARY_VASCULAR_MORPHOGENESIS",
    "GOBP_CARDIAC_VASCULAR_SMOOTH_MUSCLE_CELL_DIFFERENTIATION",
    "GOBP_ADULT_HEART_DEVELOPMENT"
  ))
selected_gene_sets
geneSets<-lapply(unique(selected_gene_sets$gs_name),
                 function(x){selected_gene_sets$gene_symbol[selected_gene_sets$gs_name==x]})
geneSets <- list(inflammatory_response = geneSets[[1]],
                 acute_inflammatory_response = geneSets[[2]],
                 chronic_inflammatory_response = geneSets[[3]],
                 lipid_biosynthetic_process = geneSets[[4]],
                 lipid_catabolic_process= geneSets[[5]],
                 lipid_oxidation= geneSets[[6]],
                 oxidative_phosphorylation= geneSets[[7]],
                 cell_proliferation_involved_in_heart_morphogenesis= geneSets[[8]],
                 cell_migration_involved_in_heart_development= geneSets[[9]],
                 sprouting_angiogenesis= geneSets[[10]],
                 apoptotic_signaling_pathway= geneSets[[11]],
                 response_to_fluid_shear_stress= geneSets[[12]],
                 cell_cell_adhesion= geneSets[[13]],
                 cell_adhesion_involved_in_heart_morphogenesis= geneSets[[14]],
                 cell_surface_receptor_signaling_pathway_involved_in_heart_development= geneSets[[15]],
                 angiogenesis_involved_in_coronary_vascular_morphogenesis= geneSets[[16]],
                 cardiac_vascular_smooth_muscle_cell_differentiation= geneSets[[17]],
                 adult_heart_development= geneSets[[18]]
)
a <- list(
  m6A = m6a_features
)
list <- c(a,geneSets)
library(GSVA)
library(corrplot)
gsva_mat <- gsva(expr=rt,
                 gset.idx.list=list,
                 kcdf="Gaussian" ,#"Gaussian" for logCPM,logRPKM,logTPM, "Poisson" for counts
                 verbose=T,
                 parallel.sz = 1)#调用所有核
saveRDS(gsva_mat,"gsva.rds")
gsva_mat <- readRDS("gsva.rds")
gsva_mat[1:5,1:5]
gsva_mat <- t(gsva_mat)
gsva_mat[1:5,1:5]
meta <- sce@meta.data
meta$name <- rownames(meta)
meta <- meta[,c(18,19)]
newscore <- cbind(meta,gsva_mat)
cell_type <- c("Fibroblast 1","Endothelial","Macrophage","Fibromyocyte","T cell","Smooth muscle cell",
               "Pericyte 1","Pericyte 2","B cell","Plasma cell 1","Fibroblast 2","Neuron","Plasma cell 2",
               "NK cell","Mast cell")
new <- data.frame('names','variable','r','p','sig')
colnames(new) <- c('names','variable','r','p','sig')
for(i in cell_type){
  celltype <- subset(newscore,newscore$cell_type == i)
  celltype <- celltype[,-c(1,2)]
  cor <- corr.test(celltype[,2:19],celltype[,1], method = "spearman", adjust = "fdr")
  cp <- as.data.frame(cor$p)
  cr <- as.data.frame(cor$r)
  cr <- round(cr,2)
  cp$names <- rownames(cp)
  longp <- melt(cp,idvar = "names",v.names = "abd",direction = "long")
  colnames(longp) <- c('na','va','p')
  cr$names <- rownames(cr)
  longr <- melt(cr,idvar = "names",v.names = "abd",direction = "long")
  allnew <- as.data.frame(cbind(longr,longp$p))
  colnames(allnew) <- c('names','variable','r','p')
  allnew[which(allnew$p<0.001),'sig'] <- '***'
  allnew[which(allnew$p<0.01 & allnew$p>0.001),'sig'] <- '**'
  allnew[which(allnew$p<0.05 & allnew$p>0.01),'sig'] <- '*'
  allnew$variable <- i
  new <- rbind(allnew,new)
}
new <- new[-271,]
new$r <- as.numeric(new$r)
new$p <- as.numeric(new$p)
class(new$r)
write_csv(new,"cordata.csv")
setwd("~/20240103_Atherosis/v2/result/Fig2/heatmap")
new <- read_csv("cordata.csv")
new$names <- factor(new$names,levels = c("adult_heart_development",
                                         "cardiac_vascular_smooth_muscle_cell_differentiation",
                                         "angiogenesis_involved_in_coronary_vascular_morphogenesis",
                                         "cell_surface_receptor_signaling_pathway_involved_in_heart_development",
                                         "cell_adhesion_involved_in_heart_morphogenesis",
                                         "cell_cell_adhesion",
                                         "response_to_fluid_shear_stress",
                                         "apoptotic_signaling_pathway",
                                         "sprouting_angiogenesis",
                                         "cell_migration_involved_in_heart_development",
                                         "cell_proliferation_involved_in_heart_morphogenesis",
                                         "oxidative_phosphorylation",
                                         "lipid_oxidation",
                                         "lipid_catabolic_process",
                                         "lipid_biosynthetic_process",
                                         "chronic_inflammatory_response",
                                         "acute_inflammatory_response",
                                         "inflammatory_response"
))
ggplot(new, aes(variable,names,fill = r)) +
  scale_fill_gradient2(low = '#864b76', high='#e23b54',mid = 'white',
                       limit=c(-0.4,0.4),name=paste0("Correlation")) +
  geom_tile(color = "white",lwd = 0.5,linetype = 1) +
  geom_text(aes(label=r), color="black", size=2) +
  geom_text(aes(label=sig), color="black", size=2,vjust = 1.8) +
  labs(x=NULL,y=NULL) +
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(size=5,angle = -45,hjust = 1,color = "black"),
        axis.text.y = element_text(size=5,color = "black"),
        axis.ticks.y = element_blank(),
        panel.background=element_blank())

ggsave("cor6_gsva.pdf",height = 6,width = 8)

# 通路气泡图-----
library(Seurat)
library(clusterProfiler)
library(tidyverse)
library(AUCell)
library(Seurat)
library(tidyverse)
source("custom_function.R")
path <-
  "~/20240103_Atherosis/v2/result/1-dealdata/seurat_integration_anno2.rds"
seurat_obj <- readRDS(path)
dim(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- seurat_obj %>%
  cat_construct_metacells(k = 15,name = "m6a_cor")
saveRDS(seurat_obj,"k15-metacell_obj.rds")
metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj)
saveRDS(metacell_obj,"k15-metacell_obj2.rds")
metacell_obj <- readRDS("~/20240103_Atherosis/v2/result/Fig2/supply/1-m6agene/k15-metacell_obj2.rds")
m6a_features <- read_csv("m6a_genesets.csv") %>%
  filter(feature %in% rownames(metacell_obj)) %>%
  pull(feature)
geneSets <- readRDS("geneSets.rds")
s_sets <- geneSets
library(AUCell)
score_output <- file.path("aucell_metacell.rds")
cell_types <- c("B_cell","C10_S100A8","C9_APOE","Endothelial","Fibroblast_1","Fibroblast_2","Fibromyocyte",
                "Macrophage", "Neuron","Smooth_muscle_cell","T_cell", "Plasma_cell_1","Plasma_cell_2",
                "Pericyte_1","Pericyte_2")
metacell_obj$cell_type <- gsub(" ","_",metacell_obj@meta.data$cell_type)
metacell_obj$cell_type_intermediate <- metacell_obj$cell_type
table(metacell_obj$cell_type)
res <-
  map(
    cell_types,
    run_aucell,
    seurat_obj = metacell_obj,
    s_sets = s_sets,
    assay = "RNA",
    score_output = score_output
  )
names(res) <- cell_types
saveRDS(res,"~/20240103_Atherosis/v2/result/Fig2/supply/2-m6apathway/aucCelltype/res.rds")
es_matrix_list <- res
cor_output <- file.path("~/20240103_Atherosis/v2/result/Fig2/supply/2-m6apathway/")
m6a_features <- c("FTO","METTL3","METTL14","RBM15","RBM15B", "WTAP","CBLL1","ZC3H13","ALKBH5",
                  "YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","IGF2BP1", "IGF2BP2","IGF2BP3",
                  "HNRNPA2B1","HNRNPC", "FMR1","LRPPRC","ELAVL1","VIRMA")
res <- tidyr::crossing(cell_types, m6a_features) %>%
  dplyr::rename(cell_type = cell_types, feature = m6a_features) %>%
  pmap_df(
    run_cor_fun,
    seurat_obj = metacell_obj,
    es_matrix = es_matrix_list,
    av = F
  )
write_csv(res,file = file.path("cor_metacell.csv"))
res <- read_csv("~/20240103_Atherosis/v2/result/Fig2/supply/2-m6apathway/cor_metacell.csv")
filter_res  <- res  %>%
  filter(feature_x %in% c("FTO","METTL3","METTL14", "WTAP","ALKBH5",
                          "YTHDC1","YTHDF1","YTHDF2","IGF2BP1", "IGF2BP3"
  ),
  p_value <= 0.05,
  abs(estimate)>0.1,
  ) %>%
  filter(cell_type %in% c("Endothelial","Macrophage","Smooth_muscle_cell",
                          "Fibromyocyte","T_cell"))
min(filter_res$p_value)
min(filter_res$estimate)
max(filter_res$estimate)
filter_res$p_value[which(filter_res$p_value<1e-20)]=1e-20
colnames(filter_res)[4] <- "cor"
df <- filter_res
df$feature_y <- gsub("GOBP_", "",df$feature_y )
df$feature_y  <- gsub("_", " ",df$feature_y )
df$feature_y  <- tolower(df$feature_y )
df$feature_y <- factor(df$feature_y,levels = c("adult heart development",
                                               "cardiac vascular smooth muscle cell differentiation",
                                               "angiogenesis involved in coronary vascular morphogenesis",
                                               "cell surface receptor signaling pathway involved in heart development",
                                               "cell adhesion involved in heart morphogenesis",
                                               "cell cell adhesion",
                                               "response to fluid shear stress",
                                               "apoptotic signaling pathway",
                                               "sprouting angiogenesis",
                                               "cell migration involved in heart development",
                                               "cell proliferation involved in heart morphogenesis",
                                               "oxidative phosphorylation",
                                               "lipid oxidation",
                                               "lipid catabolic process",
                                               "lipid biosynthetic process",
                                               "chronic inflammatory response",
                                               "acute inflammatory response",
                                               "inflammatory response"
))
p3 <-df  %>%
  ggplot() +
  geom_point(aes(x = feature_x,
                 y = feature_y,
                 fill = cor,
                 size = -log10(p_value)),shape = 22,color = "#4d474d",stroke = 0.2)+ # 23 stroke =0.2 
  facet_grid(. ~ cell_type,
             space = "free",
             scales = "free" 
  ) +
  scale_fill_gradient2(
    low = "#864b76", mid = "white", high = "#af2934")+
  theme_bw()+
  theme(
    legend.position = "top",
    panel.spacing.x = unit(0, "pt"),
    panel.grid=element_blank(),
    axis.text.x = element_text(angle = -90,hjust = 1,vjust = 1,colour = "black",size = 6),
    axis.text.y = element_text(colour = "black",size = 5),
    strip.background = element_blank()
  )
ggsave("4.pdf",p3,height = 6,width = 10)

# 热图-----
library(Seurat)
library(clusterProfiler)
library(tidyverse)
library(AUCell)
source("custom_function.R")
source("custom_plot_function.R")
metacell_obj <- readRDS("~/20240103_Atherosis/v2/result/Fig2/supply/1-m6agene/k15-metacell_obj2.rds")
dim(metacell_obj)
table(metacell_obj$cell_type)
metacell_obj$cell_type <- gsub(" ","_",metacell_obj@meta.data$cell_type)
cell_types <- c("Endothelial","Fibromyocyte","Macrophage","Smooth_muscle_cell","T_cell")
C5_gene_sets <- msigdbr::msigdbr(species = "human",
                                 category = "C5") %>%dplyr::select(gs_name, gene_symbol)
a <- as.data.frame(unique(C5_gene_sets$gs_name))
selected_gene_sets <- C5_gene_sets %>%filter(gs_name %in% c("GOBP_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY_INVOLVED_IN_HEART_DEVELOPMENT"))
selected_gene_sets
genes <- selected_gene_sets$gene_symbol 
seurat_genes <- rownames(metacell_obj@assays$RNA@data)
genes_present <- genes %in% seurat_genes
table(genes_present)
genes_not_present <- genes[!(genes %in% seurat_genes)]
genes_not_present
a <- unique(subset(genes,!genes %in% genes_not_present))
ASgenes <- read_delim("~/20240103_Atherosis/v2/result/Fig2/supply/3-m6Atermgene/cell_cell_adhesion/genecard_AS_relatedGenes.txt")
ASgenes <- colnames(ASgenes)
b <- subset(a, a %in% ASgenes)
b
selected_genes <- c("ACVR1", "TGFB1","BMP2","CTNNB1","DLL4","JAG1","NOTCH1","NOTCH2")#Tight Junction Proteins（紧密连接蛋白）
m6a_features <- read_csv("/home/pingxr/scrna/m6A/20220415_m6A-brain/data/m6a_genesets.csv") %>%
  filter(feature %in% rownames(metacell_obj)) %>%
  pull(feature)
df <- data.frame(v1 = "feature_x", v2 = "feature_y",v3 = "p_value",v4 = "estimate",v5 = "num",v6 = "cell_type")
colnames(df) <- c("feature_x","feature_y","p_value","estimate","num","cell_type")
df <- df[-1,]
for ( i in cell_types){
  genes_cor_res <- tidyr::crossing(i, m6a_features) %>%
    dplyr::rename(cell_type = i, feature_x = m6a_features) %>%
    pmap_df(
      run_cor,
      seurat_obj = metacell_obj,
      feature_y = selected_genes
    )
  head(genes_cor_res)
  write.csv(genes_cor_res,paste0(i,"_genescor.csv"))
  df <- rbind(df,genes_cor_res)
}
write.csv(df,"all_genescor.csv")
df <- read_csv("all_genescor.csv")
df$genename <- paste0(df$cell_type,"-",df$feature_y)
df$Gene.type <- "Cell_surface_receptor_signaling_pathway_involved_in_heart_development"
df1 <- df[,c(7,8,9)] 
df1 <- subset(df1, !duplicated(df1))
df2 <- df[,c(2,5,8)]
df2 <-df2 %>% 
  pivot_wider(names_from = genename, values_from = estimate)
a <- df2$feature_x
df2 <- df2[,-1]
rownames(df2) <- a 
metacell_obj <- readRDS("~/20240103_Atherosis/v2/result/Fig2/supply/1-m6agene/k15-metacell_obj2.rds")
m6a_features <- read_csv("m6a_genesets.csv") %>%
  filter(feature %in% rownames(metacell_obj))  
matrix_2 <-data.frame(scale(df2,center = T)) #中心化
annotation_col<- data.frame(Celltype = df1$cell_type,
                            Pathway_Type = df1$Gene.type,# 构建行注释信息
                            row.names = df1$genename)
row.names(annotation_col) <- colnames(df2)
type_color <- c("#F89C74","#75ACC3","#66C5CC","#80BA5A","#FE88B1")
names(type_color) <- c("Endothelial","Fibromyocyte","Macrophage","Smooth_muscle_cell","T_cell") 
Genetype_color <- c("#C2C0A6")
names(Genetype_color) <- c("Cell_surface_receptor_signaling_pathway_involved_in_heart_development") 
annotation_row <- data.frame(m6A_type = m6a_features$group,
                             row.names = m6a_features$feature)
row.names(annotation_row) <- row.names(df2)
m6Atype_color <- c("#A8B545","#6A8C69","#53736A")
names(m6Atype_color) <- c("Eraser","Reader","Writer") 
ann_colors <- list(Celltype=type_color, Pathway_Type = Genetype_color,m6A_type=m6Atype_color) 
matrix_2[is.na(matrix_2)] <- 0
colnames(matrix_2) <- gsub("\\.", "-", colnames(matrix_2) )
bk <- c(seq(-4,-0.1,by=0.01),seq(0,4,by=0.01))
pdf("1.pdf",height = 5,width = 10)
p <- pheatmap(matrix_2,
              scale="row",
              color =  colorRampPalette(c("#348888", "white", "#F24405"))(100),
              annotation_col =  annotation_col,
              annotation_row = annotation_row,
              annotation_colors = ann_colors,
              fontsize_col = 7,
              cluster_rows = T,
              treeheight_row = 0.5,
              cluster_cols = F,
              show_rownames =T, 
              show_colnames = T,
              fontsize = 7,
              cellwidth=7,
              cellheight=7, 
              main = "Correlation of the m6A gene with \nCell_surface_receptor_signaling_pathway_involved_in_heart_development") # main参数添加主标题
dev.off()
print(p)

