#### Construct_metacells ----
cat_construct_metacells <- function(seurat_obj, k, name, min_cells = 50, max_shared = 10, assay = "RNA") {
  if ("cell_type" %in% colnames(seurat_obj@meta.data) &
      "donor" %in% colnames(seurat_obj@meta.data)) {
    set.seed(717)
    seurat_obj <- hdWGCNA::SetupForWGCNA(
      seurat_obj,
      gene_select = "fraction",
      fraction = 0.05,
      wgcna_name = name
    )
    seurat_obj <- hdWGCNA::MetacellsByGroups(
      seurat_obj = seurat_obj,
      group.by = c("cell_type", "donor"),
      k = k,
      min_cells = min_cells,
      max_shared = max_shared,
      mode = "sum",
      assay = assay,
      ident.group = "cell_type"
    )
    seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)
    return(seurat_obj)
  }
  stop("The column name of cell_type or donor does not exist in the meta.data of your Seurat object, please add it!")
}
#### Correlation ----
compute_correlation <-
  function(adata,
           feature_x,
           feature_y,
           method = "pearson",
           orientation = "column") {
    if (orientation == "column") {
      cor.test.res <-
        cor.test(adata[, feature_x],
                 adata[, feature_y],
                 method = method)
      p.value <- cor.test.res$p.value
      estimate <- cor.test.res$estimate
      return(
        data.frame(
          feature_x = feature_x,
          feature_y = feature_y,
          p_value = p.value,
          estimate = estimate,
          num = nrow(adata),
          method = method
        )
      )
    }
    if (orientation == "row") {
      cor.test.res <-
        cor.test(as.numeric(adata[feature_x,]),
                 as.numeric(adata[feature_y,]),
                 method = method)
      p.value <- cor.test.res$p.value
      estimate <- cor.test.res$estimate
      return(
        data.frame(
          feature_x = feature_x,
          feature_y = feature_y,
          p_value = p.value,
          estimate = estimate,
          num = ncol(adata),
          method = method
        )
      )
    }
  }

cat_correlation <- function(seurat_obj,
                            feature_x,
                            feature_y = NULL,
                            cell_types = NULL,
                            score_method = "aucell",
                            correlation_method = "pearson",
                            species = "human",
                            ncores = 10,
                            outdir = NULL) {
  if (!is.null(cell_types) & length(cell_types) == 1) {
    seurat_obj <- subset(seurat_obj, cell_type %in% cell_types)
  }
  data_matrix <-
    as.data.frame(Seurat::GetAssayData(seurat_obj, slot = "data", assay = "RNA"))
  if (all(feature_y %in% rownames(seurat_obj)) &
      length(feature_y) >= 10) {
    temp_adata <- tidyr::crossing(feature_x, feature_y)
    x <- temp_adata$feature_x
    y <- temp_adata$feature_y
    library(foreach)
    library(doSNOW)
    library(tcltk)
    cls <- makeSOCKcluster(ncores)
    registerDoSNOW(cls)
    pb <- txtProgressBar(max = length(x),
                         style = 3,
                         char = "*")
    progress <- function(n) {
      setTxtProgressBar(pb, n)
    }
    opts <- list(progress = progress)
    genes_cor_res <- foreach(
      x = x,
      y = y,
      .combine = "rbind",
      .export = "compute_correlation",
      .packages = c("stats", "base", "sp"),
      .options.snow = opts
    ) %dopar% compute_correlation(
      adata = data_matrix,
      feature_x = x,
      feature_y = y,
      method = correlation_method,
      orientation = "row"
    )
    close(pb)
    stopCluster(cls)
    return(genes_cor_res)
  } else if (length(feature_y) == 1) {
    if (feature_y != "proportions") {
      scores_matrix <-
        cat_score(
          seurat_obj,
          category = feature_y,
          method = score_method,
          species = species,
          ncores = ncores
        )
      if (!is.null(outdir)) {
        colnames(x = scores_matrix) <-
          gsub(
            pattern = "\\.",
            replacement = "#",
            x = colnames(x = scores_matrix)
          )
        saveRDS(scores_matrix, file.path(outdir, paste0(feature_y, ".scores.rds")))
      }
      data_matrix <- data_matrix[feature_x,]
      if (any(grepl(pattern = "#", x = colnames(data_matrix))) ||
          any(grepl(pattern = "#", x = rownames(scores_matrix)))) {
        warning(
          "Cell id cannot have hashtag ('#'), replacing with dot ('.')",
          call. = FALSE,
          immediate. = TRUE
        )
        colnames(x = data_matrix) <-
          gsub(
            pattern = "#",
            replacement = "\\.",
            x = colnames(x = data_matrix)
          )
        colnames(x = scores_matrix) <- gsub(
          pattern = "#",
          replacement = "\\.",
          x = colnames(x = scores_matrix)
        )
      }
      if (identical(colnames(data_matrix), colnames(scores_matrix))) {
        data_matrix <- rbind(data_matrix, scores_matrix)
        feature_y <- rownames(scores_matrix)
        temp_adata <- tidyr::crossing(feature_x, feature_y)
        x <- temp_adata$feature_x
        y <- temp_adata$feature_y
        library(foreach)
        library(doSNOW)
        library(tcltk)
        cls <- makeSOCKcluster(ncores)
        registerDoSNOW(cls)
        pb <- txtProgressBar(max = length(x),
                             style = 3,
                             char = "*")
        progress <- function(n) {
          setTxtProgressBar(pb, n)
        }
        opts <- list(progress = progress)
        scores_cor_res <- foreach(
          x = x,
          y = y,
          .combine = "rbind",
          .export = "compute_correlation",
          .packages = c("stats", "base", "sp"),
          .options.snow = opts
        ) %dopar% compute_correlation(
          adata = data_matrix,
          feature_x = x,
          feature_y = y,
          method = correlation_method,
          orientation = "row"
        )
        close(pb)
        stopCluster(cls)
        return(scores_cor_res)
      }
    }
    if (feature_y == "proportions" & !is.null(cell_types)) {
      x_cell <- cell_types[1]
      y_cell <- cell_types[-1]
      
      x_data <-
        AverageExpression(
          subset(seurat_obj, cell_type == x_cell),
          assays = "RNA",
          slot = "data",
          group.by = "donor",
          features = feature_x
        )[["RNA"]]
      x_data <- as.data.frame(x_data)
      sub_seurat_obj <- subset(seurat_obj, cell_type %in% y_cell)
      
      if (length(y_cell) == 1) {
        y_data <- as.data.frame.matrix(prop.table(
          table(sub_seurat_obj$cell_type, sub_seurat_obj$donor),
          margin = 1
        ))
      } else {
        y_data <- as.data.frame.matrix(prop.table(
          table(sub_seurat_obj$cell_type, sub_seurat_obj$donor),
          margin = 2
        ))
      }
      
      if (identical(colnames(x_data), colnames(y_data))) {
        adata <- rbind(x_data, y_data)
        if (!is.null(outdir)) {
          saveRDS(adata, file.path(outdir, paste0(feature_y, ".data.rds")))
        }
        feature_y <- rownames(y_data)
        temp_adata <- tidyr::crossing(feature_x, feature_y)
        x <- temp_adata$feature_x
        y <- temp_adata$feature_y
        
        cell_prop_cor_res <- foreach(
          x = x,
          y = y,
          .combine = "rbind",
          .export = "compute_correlation",
          .packages = c("stats", "base", "sp")
        ) %do% compute_correlation(
          adata = adata,
          feature_x = x,
          feature_y = y,
          method = correlation_method,
          orientation = "row"
        )
        cell_prop_cor_res$feature_x <- paste0(x_cell, "_", x)
        return(cell_prop_cor_res)
      }
    }
  }
}

#### Score ----
cat_score <-
  function(object,
           category,
           method = "aucell",
           species = "human",
           ncores = 10) {
    if (class(object)[1] == "Seurat") {
      data_matrix <-
        as.matrix(Seurat::GetAssayData(object, assay = "RNA", slot = "data"))
    } else {
      data_matrix <- as.matrix(object)
    }
    
    collections <- msigdbr::msigdbr_collections()
    gs_cat <- names(table(collections$gs_cat))
    gs_subcat <- names(table(collections$gs_subcat))[-1]
    if (all(category %in% gs_cat)) {
      gene_sets <- map_dfr(category, function(x) {
        msigdbr::msigdbr(species = species,
                         category = x) %>%
          dplyr::select(gs_name, gene_symbol)
      })
    } else if (all(category %in% gs_subcat)) {
      gene_sets <- map_dfr(category, function(x) {
        msigdbr::msigdbr(species = species,
                         subcategory = x) %>%
          dplyr::select(gs_name, gene_symbol)
      })
    } else {
      stop("The category name is wrong!")
    }
    gene_sets <- split(gene_sets$gene_symbol, gene_sets$gs_name)
    
    if (method == "aucell") {
      cells_rankings <-
        AUCell::AUCell_buildRankings(data_matrix,
                                     nCores = ncores,
                                     plotStats = F)
      cells_AUC <- AUCell::AUCell_calcAUC(gene_sets, cells_rankings)
      score_adata <- data.frame(AUCell::getAUC(cells_AUC))
    }
    
    if (method == "ssgsea") {
      score_adata <-
        GSVA::gsva(
          data_matrix,
          gene_sets,
          method = c("ssgsea"),
          kcdf = c("Gaussian"),
          parallel.sz = ncores
        )
      score_adata <- data.frame(score_adata)
    }
    
    if (method == "gsva") {
      score_adata <-
        GSVA::gsva(
          data_matrix,
          gene_sets,
          method = c("gsva"),
          kcdf = c("Gaussian"),
          parallel.sz = ncores
        )
      score_adata <- data.frame(score_adata)
    }
    return(score_adata)
  }

#### Differentially expressed gene analysis ----
cat_deg <- function(
  seurat_obj,
  group_by,
  control_group,
  treatment_group,
  cell_types,
  test_use = "",
  ncores = 10,
) {
  FindMarkers(seurat_obj, )
}
