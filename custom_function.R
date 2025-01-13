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

#### Visualization ----
# 中国传统色配色
cc <- c(
  "#a91e32",
  "#e94929",
  "#f4abac",
  "#eab946",
  "#c36625",
  "#f0c800",
  "#a8cd34",
  "#62be9d",
  "#bee0d0",
  "#004ea2",
  "#le7fa0",
  "#c3e0e7",
  "#672e1d",
  "#874154",
  "#e0c4ce",
  "#ad5942",
  "#a58262",
  "#d5b210",
  "#141722",
  "#434c50",
  "#f0f0f4"
)
catlabel <- function(adata, flip = F) {
  adata <-
    p <- adata %>%
    rename_at(1, ~ "x") %>%
    arrange_at(2) %>%
    mutate(x = as_factor(x)) %>%
    ggplot() +
    geom_tile(aes_string(
      x = "x",
      y = 1,
      fill = colnames(adata)[2]
    )) +
    theme_void() +
    theme(
      strip.text = element_text(size = 8, colour = "black"),
      legend.title = element_text(size = 8, colour = "black"),
      legend.text = element_text(size = 8, colour = "black"),
      legend.key.height = unit(8, "pt"),
      legend.key.width = unit(8, "pt")
    )
  if (flip) {
    p <- p + coord_flip()
  }
  return(p)
}

catdotplot <- function(adata,
                       x,
                       y,
                       color,
                       size,
                       title = NULL,
                       annotation_row = NULL,
                       annotation_col = NULL,
                       cluster_row = FALSE,
                       cluster_col = FALSE,
                       dot_scale = 6) {
  adata <- adata %>%
    dplyr::rename(x = {
      {
        x
      }
    }, y = {
      {
        y
      }
    })
  if (!is.null(annotation_col)) {
    x_order <- annotation_col %>% pull(1)
    adata <- adata %>%
      mutate(x = as_factor(x)) %>%
      mutate(x = fct_relevel(x, x_order))
  }
  if (!is.null(annotation_row)) {
    y_order <- annotation_row %>% pull(1)
    adata <- adata %>%
      mutate(y = as_factor(y)) %>%
      mutate(y = fct_relevel(y, y_order))
  }
  p <- adata %>%
    ggplot(aes(
      x = x,
      y = y,
      size = {
        {
          size
        }
      },
      color = {
        {
          color
        }
      }
    )) +
    geom_point() +
    theme_cat() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        colour = "black",
        size = 8,
        # face = "italic"
      ),
      axis.text.y = element_text(colour = "black",
                                 # face = "italic",
                                 size = 8),
      axis.title = element_blank(),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.size = unit(8, "pt"),
      legend.key = element_blank(),
      legend.margin = margin(r = 0, unit = "pt"),
      plot.title = element_text(
        size = 8,
        colour = "black",
        face = "plain",
        hjust = 0.5
      )
    ) +
    scale_color_gradient2(low = "navy", mid = "white", high = "firebrick3") +
    scale_radius(range = c(0, dot_scale))
  if (!is.null(title)) {
    p <- p + labs(title = title)
  }
  if (!is.null(annotation_col)) {
    pc <- catlabel(annotation_col) +
      scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::Set3"))
    
    p <- p %>%
      aplot::insert_top(pc, height = 0.02)
  }
  if (!is.null(annotation_row)) {
    pr <- catlabel(annotation_row, flip = T)
    # scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::Paired"))
    p <- p %>%
      aplot::insert_right(pr, width = 0.02)
  }
  if (cluster_row) {
    cluster_row <- adata %>%
      ungroup() %>%
      dplyr::select(x, y, {
        {
          color
        }
      }) %>%
      pivot_wider(names_from = x, values_from = {
        {
          color
        }
      }) %>%
      replace(is.na(.), 0) %>%
      column_to_rownames("y") %>%
      dist() %>%
      hclust() %>%
      ggtree::ggtree()
    p <- p %>%
      aplot::insert_left(cluster_row, width = 0)
  }
  if (cluster_col) {
    cluster_col <- adata %>%
      ungroup() %>%
      dplyr::select(x, y, {
        {
          color
        }
      }) %>%
      pivot_wider(names_from = y, values_from = {
        {
          color
        }
      }) %>%
      replace(is.na(.), 0) %>%
      column_to_rownames("x") %>%
      dist() %>%
      hclust() %>%
      ggtree::ggtree() +
      ggtree::layout_dendrogram()
    p <- p %>%
      aplot::insert_top(cluster_col, height = 0)
  }
  return(p)
}

catscatter <-
  function(adata,
           x,
           y,
           color = NULL,
           size = NULL,
           method = NULL,
           cutoff_x = 0,
           cutoff_y = 0) {
    sub_adata <-
      adata %>% dplyr::select({
        {
          x
        }
      }, {
        {
          y
        }
      }, {
        {
          color
        }
      }, {
        {
          size
        }
      })
    if (!is.null(cutoff_x)) {
      sub_adata <- sub_adata %>%
        dplyr::filter({
          {
            x
          }
        } > cutoff_x)
    }
    if (!is.null(cutoff_y)) {
      sub_adata <- sub_adata %>%
        dplyr::filter({
          {
            y
          }
        } > cutoff_y)
    }
    p <- sub_adata %>%
      ggplot2::ggplot(ggplot2::aes(
        x = {
          {
            x
          }
        },
        y = {
          {
            y
          }
        },
        color = {
          {
            color
          }
        },
        size = {
          {
            size
          }
        }
      )) +
      ggrastr::rasterise(ggplot2::geom_point(size = 0.2), dpi = 600) +
      theme_cat()
    if (!is.null(method)) {
      cor_res <-
        dplyr::pull(sub_adata, {
          {
            x
          }
        }) %>% cor.test(dplyr::pull(sub_adata, {
          {
            y
          }
        }), method = method)
      p_value <- cor_res$p.value
      estimate <- cor_res$estimate
      p <- p +
        ggplot2::geom_smooth(
          method = "lm",
          formula = y ~ x,
          color = "#6b76ae",
          fill = "#cbc9e2",
          size = 0.5
        ) +
        ggplot2::labs(title = paste0(
          "r = ",
          round(estimate, 2),
          ", p = ",
          format(p_value, digits = 2, scientific = T)
        ))
    }
    return(p)
  }

theme_cat <- function(...) {
  theme(
    # Axis
    axis.title = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    axis.text = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(size = 0.5, colour = "black"),
    # Panel
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      size = 0.5,
      fill = NA,
      colour = "black"
    ),
    # Legend
    legend.background = element_blank(),
    legend.title = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    legend.text = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    legend.key = element_blank(),
    legend.key.height = unit(8, "pt"),
    legend.key.width = unit(8, "pt"),
    # Plot
    plot.background = element_blank(),
    plot.title = element_text(
      size = 8,
      hjust = 0.5,
      face = "plain",
      colour = "black"
    ),
    # plot.margin = margin(
    #   t = 0,
    #   r = 0,
    #   b = 0,
    #   l = 0
    # ),
    # Facetting
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    ...
  )
}

#### Compute ----
# 关注分子与各个细胞类型功能术语的相关性
feature_cor_function <- function(seurat_obj,
                                 cell_type,
                                 features,
                                 top_n_term = 50) {
  ident.1 <- cell_type
  print(ident.1)
  Idents(seurat_obj) <- "cell_type"
  markers <- FindMarkers(
    seurat_obj,
    ident.1 = ident.1,
    only.pos = T,
    assay = "RNA",
    slot = "data",
    logfc.threshold = 0.5
  ) %>%
    filter(p_val_adj <= 0.05)
  dim(markers)
  ego <- clusterProfiler::enrichGO(
    rownames(markers),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
  )
  res <- ego@result %>%
    filter(p.adjust <= 0.05) %>%
    top_n(top_n_term, wt = -p.adjust)
  gene_list <- str_split(res$geneID, "/")
  names(gene_list) <-
    res$Description %>%
    str_replace_all(" ", "_") %>%
    str_replace_all("'", "") %>%
    str_replace_all("-", "")
  
  sub_seurat_obj <- seurat_obj %>%
    subset(idents = ident.1)
  sub_seurat_obj <-
    AddModuleScore(
      sub_seurat_obj,
      assay = "RNA",
      features = gene_list,
      seed = 717,
      name = str_c(ident.1, "__", "gobp_", names(gene_list), "__")
    )
  expr <-
    sub_seurat_obj %>% GetAssayData(slot = "data", assay = "RNA")
  
  sub_expr <- expr[pull(m6a_genesets, feature),] %>%
    as.matrix() %>%
    t() %>%
    as.data.frame()
  sub_seurat_obj <-
    sub_seurat_obj %>% AddMetaData(sub_expr)
  
  meta_data <- sub_seurat_obj@meta.data
  gobp_features <-
    grep(
      colnames(meta_data),
      pattern = str_c("^", ident.1, "__", "gobp_"),
      value = T
    )
  
  res <- crossing(features, gobp_features) %>%
    dplyr::rename(feature.x = 1, feature.y = gobp_features) %>%
    pmap_df(calculate_correlation,
            adata = meta_data,
            cutoff.x = 0) %>%
    mutate(cell_type = ident.1) %>%
    mutate(order = str_split_fixed(feature_y, "__", 3)[, 3]) %>%
    mutate(cell_type_function = str_split_fixed(feature_y, "__", 3)[, 1]) %>%
    mutate(feature_y = str_split_fixed(feature_y, "__", 3)[, 2])
  return(res)
}
# tibble_to_list
tibble_to_list <- function(tibble, x, y) {
  init_list <- list()
  seletcted_tibble <- tibble %>%
    dplyr::rename(x = {
      {
        x
      }
    }, y = {
      {
        y
      }
    }) %>%
    dplyr::select(x, y)
  
  x_names <- names(table(seletcted_tibble$x))
  x_name <- x_names[1]
  get_one_list <- function(seletcted_tibble, x_name) {
    seletcted_tibble %>%
      filter(x == x_name) %>%
      pull(y)
  }
  
  final_list <- map(x_names, get_one_list,
                    seletcted_tibble = seletcted_tibble)
  names(final_list) <- x_names
  return(final_list)
}


calculate_correlation <-
  function(adata,
           feature_x,
           feature_y,
           cutoff_x = -1,
           cutoff_y = -1,
           method = "pearson") {
    if (!is.null(cutoff_x)) {
      adata <-
        adata[adata[, feature_x] > cutoff_x,]
    }
    if (!is.null(cutoff_y)) {
      adata <-
        adata[adata[, feature_y] > cutoff_y,]
    }
    
    
    if (!is.null(nrow(adata))) {
      if (nrow(adata) >= 3) {
        cor.test.res <-
          cor.test(adata[, feature_x],
                   adata[, feature_y],
                   method = method)
        p.value <- cor.test.res$p.value
        estimate <- cor.test.res$estimate
        return(
          tibble(
            feature_x = feature_x,
            feature_y = feature_y,
            p_value = p.value,
            estimate = estimate,
            num = nrow(adata)
          )
        )
      }
    }
  }

run_cor <-
  function(seurat_obj,
           cell_type,
           group = NULL,
           feature_x,
           feature_y,
           av = F) {
    s_cell_type <- cell_type
    s_group <- group
    print(s_cell_type)
    print(feature_x)
    if (!is.null(group)) {
      seurat_obj <- subset(seurat_obj, subset = group == s_group)
    }
    if (!is.null(s_cell_type)) {
      seurat_obj <- subset(seurat_obj, subset = cell_type == s_cell_type)
    }
    if (av) {
      expr <-
        AverageExpression(seurat_obj, assays = "RNA", group.by = "donor")[["RNA"]] %>% t()
    }
    if (!av) {
      expr <- GetAssayData(seurat_obj,
                           assay = "RNA",
                           slot = "data") %>%
        as.matrix() %>%
        t()
    }
    # features <- setdiff(colnames(expr), feature)
    # 计算
    tictoc::tic()
    res <-
      map2_dfr(feature_x,
               feature_y,
               calculate_correlation,
               adata =
                 expr) %>%
      mutate(group = group, cell_type = s_cell_type)
    tictoc::toc()
    return(res)
  }

# 测试函数
# calculate_correlation(
#   expr,
#   feature.x = interest.gene,
#   feature.y = "TAGLN",
#   cutoff.x = 0,
#   cutoff.y = 0
# )

# 测试函数
# plot_correlation(expr,
#                  "WTAP",
#                  "TAGLN",
#                  cutoff.x = 0,
#                  cutoff.y = 0)

# 测试函数
# cell_types <- names(table(seurat.obj$cell_type))

# features <- setdiff(colnames(expr), interest.gene)
# run_cor(seurat.obj, cell_types[1], "AVM", interest.gene, "TAGLN")

run_GO <-
  function(adata,
           cell_type,
           group,
           regulate,
           estimate = 0.1,
           p_value = 0.05,
           top_n = 200,
           feature_col_name,
           estimate_col_name) {
    p <- p_value
    est <- estimate
    g <- group
    cp <- cell_type
    print(cp)
    df <- df %>%
      dplyr::rename(estimate = {
        {
          estimate_col_name
        }
      })
    if (regulate == "pos") {
      genes <- df %>%
        filter(p_value <= p,
               estimate >= est,
               group == g,
               cell_type == cp) %>%
        arrange(desc(estimate)) %>%
        top_n(top_n, wt = estimate) %>%
        pull({
          {
            feature_col_name
          }
        })
    }
    if (regulate == "neg") {
      genes <- df %>%
        filter(p_value <= p,
               estimate <= -est,
               group == g,
               cell_type == cp) %>%
        arrange(estimate) %>%
        top_n(top_n, wt = -estimate) %>%
        pull({
          {
            feature_col_name
          }
        })
    }
    print(length(genes))
    ego <- clusterProfiler::enrichGO(
      genes,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = "BP",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
    )
    ego_res <- ego@result %>%
      mutate(group = group,
             regulate = regulate,
             cell_type = cp) %>%
      as_tibble()
    return(ego_res)
  }


run_cor_fun <- function(seurat_obj,
                        es_matrix,
                        feature,
                        av = F,
                        group = NULL,
                        cell_type = NULL) {
  s_cell_type <- cell_type
  s_group <- group
  # seurat_obj <- seurat.obj
  # feature <- interest.gene
  print(s_cell_type)
  print(feature)
  if (!is.null(group)) {
    seurat_obj <- subset(seurat_obj, subset = group == s_group)
  }
  if (!is.null(cell_type)) {
    seurat_obj <- subset(seurat_obj, subset = cell_type == s_cell_type)
  }
  if (!av) {
    expr <-
      GetAssayData(seurat_obj,
                   assay = "RNA",
                   slot = "data") %>%
      as.matrix() %>%
      t()
    if (is.list(es_matrix)) {
      t_es_matrix <- es_matrix[[s_cell_type]][, rownames(expr)] %>% t()
    }
    if (!is.list(es_matrix)) {
      t_es_matrix <- es_matrix[, rownames(expr)] %>% t()
    }
    features <- colnames(t_es_matrix)
    es_matrix <- cbind(t_es_matrix, expr)[, c(features, feature)]
  }
  if (av) {
    seurat_obj[["score"]] <-
      CreateAssayObject(counts = es_matrix[[s_cell_type]])
    expr1 <-
      AverageExpression(seurat_obj,
                        assays = "RNA",
                        slot = "data",
                        group.by = "donor")[["RNA"]] %>% t()
    expr2 <-
      AverageExpression(seurat_obj,
                        assays = "score",
                        slot = "counts",
                        group.by = "donor")[["score"]] %>% t()
    colnames(expr2) <- colnames(expr2) %>% str_replace_all("-", "_")
    features <- colnames(expr2)
    es_matrix <- cbind(expr1, expr2)[, c(features, feature)]
  }
  # 计算
  res <-
    map2_dfr(feature,
             features,
             calculate_correlation,
             adata =
               es_matrix)
  res <- res %>% mutate(group = s_group, cell_type = s_cell_type)
  return(res)
}

# 计算每个donor中感兴趣基因阳性细胞的每个基因平均表达量
average_expr <- function(seurat_obj, cell_type, feature) {
  s_cell_type <- cell_type
  res <- AverageExpression(
    subset(seurat_obj, cell_type == s_cell_type),
    # 计算每个donor中基因平均表达量
    assay = "RNA",
    slot = "data",
    group.by = "donor"
  )[["RNA"]] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("donor") %>%
    dplyr::select(donor, feature) %>%
    mutate(cell_type = s_cell_type)
  return(res)
}

run_prop_cor <- function(df, cell_type, feature) {
  s_cell_type <- cell_type
  df %>%
    filter(cell_type == s_cell_type) %>%
    as.data.frame() %>%
    calculate_correlation(feature, "prop", method = "pearson") %>%
    mutate(feature_y = s_cell_type)
}

run_deg <- function(file_path,
                    prefix,
                    feature = NULL,
                    group = NULL,
                    treatment_group = NULL,
                    subset_col_name = NULL,
                    subset_group = NULL,
                    test.use,
                    logfc.threshold,
                    min.pct,
                    save_results = F) {
  print(feature)
  seurat_obj <- readRDS(file_path)
  seurat_obj
  if (!is.null(subset_col_name)) {
    seurat_obj <-
      AddMetaData(seurat_obj,
                  metadata = seurat_obj[[subset_col_name]],
                  col.name = "tmp_group")
    seurat_obj <-
      subset(seurat_obj, subset = tmp_group %in% subset_group)
  }
  # 细胞类型
  cell_types <- names(table(seurat_obj$cell_type))
  # 分组
  if (!is.null(feature)) {
    pn <-
      as.data.frame(t(ifelse(
        as.matrix(seurat_obj[feature,]@assays$SCT@data) > 0,
        "postive",
        "negative"
      )))
    seurat_obj <-
      AddMetaData(seurat_obj,
                  metadata = pn,
                  col.name = "group")
    treatment_group <- "postive"
  }
  if (!is.null(group)) {
    seurat_obj <-
      AddMetaData(seurat_obj,
                  metadata = seurat_obj[[group]],
                  col.name = "group")
  }
  groups <- names(table(seurat_obj$group))
  groups
  seurat_obj$cell_type_group <-
    str_c(seurat_obj$cell_type,
          seurat_obj$group,
          sep = "_")
  # 一组中少于3个细胞的细胞类型将被剔除
  rm_cell_type <- table(seurat_obj$cell_type, seurat_obj$group) %>%
    as.data.frame.array() %>%
    rownames_to_column("cell_type") %>%
    pivot_longer(-"cell_type") %>%
    filter(value < 3) %>%
    pull(cell_type)
  rm_cell_type
  sub_cell_types <- cell_types[!(cell_types %in% rm_cell_type)]
  
  Idents(seurat_obj) <- seurat_obj$cell_type_group
  
  FindDEG <- function(seurat_obj,
                      cell_type,
                      g1,
                      g2,
                      test.use = "wilcox",
                      logfc.threshold = logfc.threshold,
                      min.pct = min.pct) {
    print(cell_type)
    deg <- FindMarkers(
      seurat_obj,
      ident.1 = str_c(cell_type, g1, sep = "_"),
      ident.2 = str_c(cell_type, g2, sep = "_"),
      test.use = test.use,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct
    ) %>%
      rownames_to_column("feature") %>%
      mutate(cell_type = cell_type)
    return(deg)
  }
  # 测试函数
  # cell_type <- cell_types[1]
  # FindDEG(seurat_obj, cell_types[2], g1 = "AVM", g2 = "CTRL")
  g1 <- treatment_group
  g2 <- groups[!(groups %in% treatment_group)]
  deg <- map_df(
    sub_cell_types,
    FindDEG,
    seurat_obj = seurat_obj,
    g1 = g1,
    g2 = g2,
    test.use = test.use,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct
  )
  print(deg)
  # if (!is.null(group)) {
  #   deg <- deg %>% mutate(group=group)
  # }
  # if (!is.null(feature)) {
  #   deg <- deg %>% mutate(feature=feature)
  # }
  if (save_results) {
    if (!dir.exists("./output/deg/")) {
      dir.create("./output/deg/")
    }
    if (!is.null(feature)) {
      save_path <- str_c("./output/deg/",
                         prefix,
                         ".",
                         feature,
                         ".",
                         g1,
                         "_vs_",
                         g2,
                         ".",
                         test.use,
                         ".csv")
    }
    if (!is.null(group)) {
      save_path <- str_c("./output/deg/",
                         prefix,
                         ".",
                         group,
                         ".",
                         g1,
                         "_vs_",
                         g2,
                         ".",
                         test.use,
                         ".csv")
    }
    if (!is.null(subset_col_name)) {
      if (length(subset_group) != 1) {
        subset_group <- str_c(subset_group, collapse = "_")
      }
      save_path <- str_c(
        "./output/deg/",
        prefix,
        ".",
        subset_col_name,
        "-",
        subset_group,
        ".",
        feature,
        ".",
        g1,
        "_vs_",
        g2,
        ".",
        test.use,
        ".csv"
      )
    }
    write_csv(deg,
              file = save_path)
  }
  if (!save_results) {
    return(deg)
  }
}


#### 打分 ----
run_aucell <-
  function(seurat_obj, cell_type, s_sets, assay = "RNA", score_output) {
    rds_output <- score_output %>%
      str_replace(".rds", str_c(".", cell_type, ".rds"))
    if (!file.exists(rds_output)) {
      selected_cell_type <- cell_type
      print(selected_cell_type)
      sub_seurat_obj <-
        subset(seurat_obj, subset = cell_type == selected_cell_type)
      exprMatrix <-
        as.matrix(GetAssayData(sub_seurat_obj, assay = assay, slot = "data"))
      cells_rankings <- AUCell_buildRankings(exprMatrix)
      
      cells_AUC <-
        AUCell_calcAUC(
          s_sets,
          cells_rankings,
          aucMaxRank = nrow(cells_rankings) *
            0.05
        )
      saveRDS(cells_AUC, rds_output)
    }
    if (file.exists(rds_output)) {
      cells_AUC <- readRDS(rds_output)
    }
    scores_matrix <- getAUC(cells_AUC)
    scores_matrix <- scores_matrix %>% as.data.frame()
    return(scores_matrix)
  }

#### hdWGCNA ----
ConstructNetwork <- function (seurat_obj, soft_power = NULL, tom_outdir = "TOM", 
                              use_metacells = TRUE, setDatExpr = TRUE, group.by = NULL, 
                              group_name = NULL, consensus = FALSE, multi.group.by = NULL, 
                              multi_groups = NULL, blocks = NULL, maxBlockSize = 30000, 
                              randomSeed = 12345, corType = "pearson", consensusQuantile = 0.3, 
                              networkType = "signed", TOMType = "signed", TOMDenom = "min", 
                              scaleTOMs = TRUE, scaleQuantile = 0.8, sampleForScaling = TRUE, 
                              sampleForScalingFactor = 1000, useDiskCache = TRUE, chunkSize = NULL, 
                              deepSplit = 4, pamStage = FALSE, detectCutHeight = 0.995, 
                              minModuleSize = 50, mergeCutHeight = 0.2, saveConsensusTOMs = TRUE, 
                              ...) 
{
  if (consensus) {
    if (!("multiExpr" %in% names(GetActiveWGCNA(seurat_obj))) | 
        setDatExpr == TRUE) {
      seurat_obj <- SetMultiExpr(seurat_obj, group_name = group_name, 
                                 group.by = group.by, multi.group.by = multi.group.by, 
                                 multi_groups = multi_groups)
    }
    multiExpr <- GetMultiExpr(seurat_obj)
    checkSets(multiExpr)
  }
  else {
    if (!("datExpr" %in% names(GetActiveWGCNA(seurat_obj))) | 
        setDatExpr == TRUE) {
      seurat_obj <- SetDatExpr(seurat_obj, group_name = group_name, 
                               group.by = group.by, use_metacells = use_metacells, 
                               return_seurat = TRUE)
    }
    datExpr <- GetDatExpr(seurat_obj)
    if (is.null(group_name)) {
      group_name <- "all"
    }
    nSets = 1
    setLabels = gsub(" ", "_", group_name)
    shortLabels = setLabels
    multiExpr <- list()
    multiExpr[[group_name]] <- list(data = datExpr)
    checkSets(multiExpr)
  }
  if (!dir.exists(tom_outdir)) {
    dir.create(tom_outdir)
  }
  if (is.null(soft_power) & !consensus) {
    soft_power <- GetPowerTable(seurat_obj) %>% subset(SFT.R.sq >= 
                                                         0.8) %>% .$Power %>% min
    cat(paste0("Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = ", 
               soft_power, "\n"))
  }
  else if (consensus) {
    power_tables <- GetPowerTable(seurat_obj) %>% dplyr::group_split(group)
    soft_power <- sapply(power_tables, function(power_table) {
      power_table %>% subset(SFT.R.sq >= 0.8) %>% .$Power %>% 
        min
    })
    cat(paste0("Soft power not provided. Automatically using the lowest power that meets 0.8 scale-free topology fit. Using soft_power = c(", 
               paste0(soft_power, collapse = ","), ")\n"))
  }
  net <- WGCNA::blockwiseConsensusModules(multiExpr, power = soft_power, 
                                          blocks = blocks, maxBlockSize = maxBlockSize, randomSeed = randomSeed, 
                                          corType = corType, consensusQuantile = consensusQuantile, 
                                          networkType = networkType, TOMType = TOMType, TOMDenom = TOMDenom, 
                                          scaleTOMs = scaleTOMs, scaleQuantile = scaleQuantile, 
                                          sampleForScaling = sampleForScaling, sampleForScalingFactor = sampleForScalingFactor, 
                                          useDiskCache = useDiskCache, chunkSize = chunkSize, 
                                          deepSplit = deepSplit, pamStage = pamStage, detectCutHeight = detectCutHeight, 
                                          minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, 
                                          saveConsensusTOMs = saveConsensusTOMs, consensusTOMFilePattern = "ConsensusTOM-block.%b.rda", 
                                          ...)
  file.rename("ConsensusTOM-block.1.rda", paste0(tom_outdir, "/", gsub(" ", 
                                                                       "_", group_name), "_ConsensusTOM-block.1.rda"))
  params <- list(power = soft_power, blocks = blocks, maxBlockSize = maxBlockSize, 
                 randomSeed = randomSeed, corType = corType, consensusQuantile = consensusQuantile, 
                 networkType = networkType, TOMType = TOMType, TOMDenom = TOMDenom, 
                 scaleTOMs = scaleTOMs, scaleQuantile = scaleQuantile, 
                 sampleForScaling = sampleForScaling, sampleForScalingFactor = sampleForScalingFactor, 
                 useDiskCache = useDiskCache, chunkSize = chunkSize, 
                 deepSplit = deepSplit, pamStage = pamStage, detectCutHeight = detectCutHeight, 
                 minModuleSize = minModuleSize, mergeCutHeight = mergeCutHeight, 
                 saveConsensusTOMs = saveConsensusTOMs)
  seurat_obj <- SetWGCNAParams(seurat_obj, params)
  net$TOMFiles <- paste0(getwd(), '/', tom_outdir, '/', group_name, '_', net$TOMFiles)
  
  seurat_obj <- SetNetworkData(seurat_obj, net)
  mods <- GetNetworkData(seurat_obj)$colors
  seurat_obj <- SetModules(seurat_obj, mod_df = data.frame(gene_name = names(mods), 
                                                           module = factor(mods, levels = unique(mods)), color = mods))
  seurat_obj
}
