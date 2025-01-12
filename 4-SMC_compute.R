# metacell-----
setwd("/home/pingxr/Atherosis_0723/20230204_Atherosis/code/SMC/")
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")

seurat_obj_path <-
  "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/seurat/seurat_integration_anno.rds"
selected_cell_type <- "SMC"  
k <- 20  
outdir <- "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/SMC/"
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "data")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
seurat_obj <- seurat_obj_path %>% readRDS()
sub_seurat_obj <- seurat_obj %>%
  subset(cell_type == selected_cell_type)
sub_seurat_obj <- sub_seurat_obj %>%
  cat_construct_metacells(k = k, name = selected_cell_type)
metacell_obj <- hdWGCNA::GetMetacellObject(sub_seurat_obj)
table(metacell_obj$donor)
sub_seurat_obj %>% saveRDS(file.path(outdir, "seurat_obj.rds"))
head(sub_seurat_obj)

# cor -----
setwd("/home/pingxr/Atherosis_0723/20230204_Atherosis/code/SMC/")
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")

seurat_obj_path <-
  "/home/pingxr/Atherosis_0723/result_human/SMC/SMC/k_20/data/seurat_obj.rds"
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2","IGF2BP3","YTHDF1","YTHDF2")
selected_cell_type <- "SMC"
k <- 20
outdir <- "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/SMC/"
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
# wgcna ----
setwd("/home/pingxr")     
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")

selected_cell_type <- "SMC"
selected_genes <- c("METTL3","METTL14","WTAP","FTO","ALKBH5","IGF2BP2",
                    "IGF2BP3","YTHDF1","YTHDF2")
k <- 20
outdir <- "./Atherosis_0723/20230204_Atherosis/result/SMC"   
outdir <-
  file.path(outdir, selected_cell_type, str_c("k_", k), "wgcna")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = T)
}
seurat_obj <- readRDS("/home/pingxr/Atherosis_0723/20230204_Atherosis/result/SMC/SMC/k_20/data/seurat_obj.rds")
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
                   soft_power = 8,
                   setDatExpr = FALSE,
                   tom_outdir = outdir,
                   nThreads = 10)
png(file.path(outdir, "plot_dendrogram.pdf"))
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
             max_genes = 200) 
saveRDS(seurat_obj, file = file.path(outdir, "wgcna.rds"))
enrichrtable <- hdWGCNA::GetEnrichrTable(seurat_obj)
enrichrtable %>%
  write_csv(file.path(outdir, "enrichr_table.csv"))

# deg ----
setwd("/home/pingxr/Atherosis_0723/20230204_Atherosis/code/SMC/")
library(Seurat)
library(tidyverse)
source("/home/tutorial/20220802_scrna-m6A/code/computing/custom_function.R")

seurat_obj_path <-
  "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/SMC/SMC/k_20/data/seurat_obj.rds"
selected_genes <- c("WTAP")
selected_cell_type <- "SMC"
k <- 20
outdir <- "/home/pingxr/Atherosis_0723/20230204_Atherosis/result/SMC/"
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
table(metacell_obj@assays$RNA@counts['WTAP', ] > median(metacell_obj@assays$RNA@counts['WTAP', ]))
write_csv(deg_metacell_obj, file.path(outdir, "deg_metacell_obj.csv"))
enriched <- deg_metacell_obj |>
  filter(group == "WTAP") |>
  cat_enrich()
saveRDS(enriched, file.path(outdir, "WTAP_enriched.rds"))
gsea_res <- deg_metacell_obj |>
  filter(group == "WTAP") |>
  cat_gsea(category = "H"
  )
saveRDS(gsea_res, file.path(outdir, "h_gsea.rds"))
gsea_res1 <- deg_metacell_obj |>
  filter(group == "WTAP") |>
  cat_gsea(category = "C2"
  )
saveRDS(gsea_res1, file.path(outdir, "C2_gsea.rds"))
gsea_res2 <- deg_metacell_obj |>
  filter(group == "WTAP") |>
  cat_gsea(category = "C5"
  )
saveRDS(gsea_res2, file.path(outdir, "C5_gsea.rds"))


