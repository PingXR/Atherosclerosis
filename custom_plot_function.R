#### 主题 ----
theme_cat <- function(...) {
  theme(
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
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      colour = "black"
    ),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      size = 0.25,
      fill = NA,
      colour = "black"
    ),
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
    plot.background = element_blank(),
    plot.title = element_text(
      size = 8,
      hjust = 0.5,
      face = "plain",
      colour = "black"
    ),
    plot.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 0
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    ...
  )
}
cc <- c(
  "#1e7fa0",
  "#e94929",
  "#62be9d",
  "#004ea2",
  "#f4abac",
  "#bee0d0",
  "#a91e32",
  "#eab946",
  "#c36625",
  "#f0c800",
  "#a8cd34",
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
catbarplot <- function(adata,
                       x,
                       y,
                       fill = NULL,
                       position = "stack",
                       width = 0.8,
                       title = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       legend = TRUE,
                       errorbar = FALSE,
                       sort = FALSE,
                       flip = FALSE,
                       group_by = NULL,
                       label_position = NULL,
                       hline = NULL,
                       jitter = FALSE) {
  x <- enquo(x)
  y <- enquo(y)
  fill <- enquo(fill)
  group_by <- enquo(group_by)

  if (sort) {
    adata <- adata %>%
      mutate({{ x }} := fct_reorder({{ x }}, {{ y }}))
  }
  if (as_label(group_by) != "NULL") {
    adata <- adata %>%
      mutate({{ x }} := fct_reorder({{ x }}, {{ group_by }}))
  }
  if (as_label(fill) %in% colnames(adata) & !errorbar) {
    p <- adata %>%
      ggplot(aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
      geom_bar(
        stat = "identity",
        width = width,
        position = position
      )
  }
  if (!(as_label(fill) %in% colnames(adata)) & !errorbar) {
    fill <- as_label(fill) %>% str_replace_all("\"", "")
    if (fill == "NULL") {
      fill <- cc[1]
    }
    p <- adata %>%
      ggplot(aes(x = {{ x }}, y = {{ y }})) +
      geom_bar(
        stat = "identity",
        width = width,
        position = position,
        fill = fill
      )
  }
  if (errorbar) {
    if (as_label({{ x }}) == as_label({{ fill }})) {
      adata <- adata %>% group_by({{ x }})
    }
    if (as_label({{ x }}) != as_label({{ fill }})) {
      adata <- adata %>% group_by({{ x }}, {{ fill }})
    }
    p <- adata %>%
      summarise(mean = mean({{ y }}), sd = sd({{ y }})) %>%
      ggplot(aes(x = {{ x }}, y = mean, fill = {{ fill }})) +
      geom_errorbar(
        mapping = aes(
          x = {{ x }},
          ymin = mean - sd,
          ymax = mean + sd
        ),
        width = 0.2,
        size = 0.5 * 0.36,
        color = "black",
        position = position_dodge(width)
      ) +
      geom_bar(
        stat = "identity",
        width = width,
        position = position
      )
  }
  if (flip) {
    p <- p + coord_flip() +
      labs(y = xlab, x = ylab, title = title)
  }
  if (!flip) {
    p <- p +
      labs(x = xlab, y = ylab, title = title)
  }
  p <- p +
    theme_cat() +
    scale_fill_manual(values = cc)
  if (!legend) {
    p <- p + guides(fill = "none")
  }
  if (!is.null(hline)) {
    p <- p + geom_hline(
      yintercept = hline,
      color = "white",
      linetype = "dashed",
      size = 0.5
    )
  }
  if (!is.null(label_position)) {
    match.arg(label_position, c("inward", "outward"))
    if (label_position == "inward") {
      adata <- adata %>%
        mutate(hjust = if_else({{ y }} > 0, "inward", "outward")) %>%
        mutate(temp_label = if_else({{ y }} > 0, str_c(" ", {{ x }}), str_c({{ x }}, " ")))
    }
    if (label_position == "outward") {
      adata <- adata %>%
        mutate(hjust = if_else({{ y }} > 0, "outward", "inward")) %>%
        mutate(temp_label = if_else({{ y }} > 0, str_c({{ x }}, " "), str_c(" ", {{ x }})))
    }
    p <- p + geom_text(
      data = adata,
      aes(
        x = {{ x }},
        y = 0,
        label = temp_label,
        hjust = hjust
      ),
      size = 6 * 0.36
    ) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  return(p)
}

#### 气泡图 ----
catlabel <- function(adata, flip = F) {
  adata <-
    p <- adata %>%
    rename_at(1, ~"x") %>%
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
      strip.text = element_text(size = 6, colour = "black"),
      legend.title = element_text(size = 6, colour = "black"),
      legend.text = element_text(size = 6, colour = "black"),
      legend.key.height = unit(6, "pt"),
      legend.key.width = unit(6, "pt")
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
  x <- enquo(x)
  y <- enquo(y)
  size <- enquo(size)
  color <- enquo(color)

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
      x = {{ x }},
      y = {{ y }},
      size = {{ size }},
      color = {{ color }}
    )) +
    geom_point() +
    theme_cat() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        colour = "black",
        size = 6,
      ),
      axis.text.y = element_text(
        colour = "black",
        size = 6
      ),
      axis.title = element_blank(),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.key.size = unit(6, "pt"),
      legend.key = element_blank(),
      legend.margin = margin(r = 0, unit = "pt"),
      plot.title = element_text(
        size = 6,
        colour = "black",
        face = "plain",
        hjust = 0.5
      )
    ) +
    scale_color_gradient2(low = "#00a6e1", mid = "white", high = "#ee6470") +
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
    p <- p %>%
      aplot::insert_right(pr, width = 0.02)
  }
  if (cluster_row) {
    cluster_row <- adata %>%
      ungroup() %>%
      select(x, y, {{ color }}) %>%
      pivot_wider(names_from = x, values_from = {{ color }}) %>%
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
      select(x, y, {{ color }}) %>%
      pivot_wider(names_from = y, values_from = {{ color }}) %>%
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
           cutoff_x = NULL,
           cutoff_y = NULL) {
    sub_adata <-
      adata %>% dplyr::select({{ x }}, {{ y }}, {{ color }}, {{ size }})
    if (!is.null(cutoff_x)) {
      sub_adata <- sub_adata %>%
        filter({{ x }} > cutoff_x)
    }
    if (!is.null(cutoff_y)) {
      sub_adata <- sub_adata %>%
        filter({{ y }} > cutoff_y)
    }
    p <- sub_adata %>%
      ggplot2::ggplot(ggplot2::aes(
        x = {{ x }},
        y = {{ y }},
        color = {{ color }},
        size = {{ size }}
      )) +
      geom_point_rast(
        size = 0.5,
        color = "grey",
        fill = "black", raster.dpi = 600
      ) +
      ggrastr::rasterise(ggplot2::geom_point(
        size = 0.5,
        color = "grey",
        fill = "black"
      ),
      dpi = 600) +
      theme_cat()
    if (!is.null(method)) {
      cor_res <-
        pull(sub_adata, {{ x }}) %>% cor.test(pull(sub_adata, {{ y }}), method = method)
      p_value <- cor_res$p.value
      estimate <- cor_res$estimate
      p <- p +
        ggplot2::geom_smooth(
          method = "lm",
          formula = y ~ x,
          color = "#6b76ae",
          fill = "#e5a323",
          size = 0.5,
          alpha = 0.2
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

#### 棒棒糖图 ----
catlollipop <- function() {

}
#### 散点图 ----
catpoint <- function(adata,
                     x,
                     y,
                     fill,
                     color,
                     label = NULL,
                     text = NULL,
                     title = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     legend = TRUE) {
  adata <- adata %>%
    dplyr::rename(
      x = {{ x }},
      y = {{ y }},
      color = {{ color }},
      label = {{ label }}
    )

  if (!is.null(text)) {
    adata <-
      adata %>% mutate(label = ifelse(label %in% text, label, ""))
  }

  p <- adata %>%
    ggplot(aes(x = x, y = y, color = color)) +
    geom_point() +
    labs(x = xlab, y = ylab, title = title) +
    theme_cat(aspect.ratio = 1)
  if (!is.null(label)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = label),
      max.overlaps = 1000000,
      color = "black",
      size = 3.5
    )
  }
  if (!legend) {
    p <- p + guides(color = "none")
  }
  return(p)
}

#### 火山图 ----
library(tidyverse)

catvolcano <-
  function(adata,
           x,
           y,
           label,
           text = NULL,
           p_value = 0.05,
           log2FC = 0.25) {
    adata <- adata %>%
      dplyr::rename(x = {{ x }}, y = {{ y }}, label = {{ label }})
    adata <- adata %>%
      mutate(change = ifelse(
        y <= p_value & abs(x) >= log2FC,
        ifelse(x >= log2FC, "Upregulated", "Downregulated"),
        "Stable"
      ))
    max_x <- ceiling(max(abs(adata$x)))
    if (!is.null(text)) {
      adata <- adata %>% mutate(label = ifelse(label %in% text, text, ""))
    }
    print(adata %>% count(change))
    p <- adata %>%
      ggplot(aes(
        x = x,
        y = -log10(y),
        colour = change
      )) +
      ggrastr::geom_point_rast(raster.dpi = 600, size = 0.5) +
      geom_vline(
        xintercept = c(-log2FC, log2FC),
        lty = 2,
        col = "black",
        lwd = 0.5
      ) +
      geom_hline(
        yintercept = -log10(p_value),
        lty = 2,
        col = "black",
        lwd = 0.5
      ) +
      labs(
        x = "log2(Fold change)",
        y = "-log10 (P value)"
      ) +
      scale_color_manual(
        values = c("#00a6e1", "lightgrey", "#ee6470"),
        labels = c(
          str_c(
            "Down (",
            adata %>% count(change) %>% filter(change == "Downregulated") %>% pull(n),
            ")"
          ),
          "Stable",
          str_c(
            "Up (",
            adata %>% count(change) %>% filter(change == "Upregulated") %>% pull(n),
            ")"
          )
        )
      ) +
      theme_cat() +
      theme(
        aspect.ratio = 1,
        legend.position = "top",
        legend.title = element_blank(),
        legend.margin = margin(b = -8)
      ) +
      xlim(c(-max_x, max_x))
    # theme(
    #   plot.title = element_blank(),
    #   legend.position = "top",
    #   legend.title = element_blank(),
    #   legend.text = element_text(size = 6, colour = "black"),
    #   legend.key = element_blank(),
    #   aspect.ratio = 1,
    #   panel.background = element_blank(),
    #   panel.grid = element_blank(),
    #   axis.text = element_text(size = 6, colour = "black"),
    #   axis.ticks = element_line(size = 0.5, colour = "black"),
    #   axis.line = element_blank(),
    #   axis.title = element_text(size = 6, colour = "black"),
    #   panel.border = element_rect(
    #     size = 0.5,
    #     colour = "black",
    #     fill = F
    #   )
    # )

    if (!is.null(text)) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = label),
        max.overlaps = 1000,
        color = "black",
        size = 3.5,
        fontface = "italic"
      )
    }

    return(p)
  }


#### GSEA ----
tableGrob2 <- function(d, p = NULL) {
  d <- d[order(rownames(d)), ]
  tp <- gridExtra::tableGrob(d,
    theme = gridExtra::ttheme_minimal(base_size = 6)
  )
  if (is.null(p)) {
    return(tp)
  }
  p_data <- ggplot_build(p)$data[[1]]
  p_data <- p_data[order(p_data[["group"]]), ]
  pcol <- unique(p_data[["colour"]])
  j <- which(tp$layout$name == "rowhead-fg")
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i + 1]][["gp"]] <- grid::gpar(col = pcol[i])
  }
  return(tp)
}
catgsea <-
  function(x,
           geneSetID,
           title = "",
           color = "#8ac45f",
           base_size = 6,
           rel_heights = c(1.6, 0.4, 1),
           subplots = 1:3,
           pvalue_table = FALSE,
           ES_geom = "line") {
    ES_geom <- match.arg(ES_geom, c("line", "dot"))
    geneList <- position <- NULL
    if (length(geneSetID) == 1) {
      gsdata <- enrichplot:::gsInfo(x, geneSetID)
    } else {
      gsdata <-
        do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
    }
    p <- ggplot(gsdata, aes_(x = ~x)) +
      xlab(NULL) +
      theme_cat() +
      theme(axis.line.x.bottom = element_blank()) +
      scale_x_continuous(expand = c(0, 0))
    if (ES_geom == "line") {
      es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description),
        size = 0.5
      )
    } else {
      es_layer <-
        geom_point(
          aes_(y = ~runningScore, color = ~Description),
          size = 0.5,
          data = subset(gsdata, position == 1)
        )
    }
    p.res <- p + es_layer + theme(
      legend.position = c(0.8, 0.8),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    )
    p.res <- p.res + ylab("Enrichment Score") + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0,
        unit = "cm"
      )
    )
    i <- 0
    for (term in unique(gsdata$Description)) {
      idx <- which(gsdata$ymin != 0 & gsdata$Description ==
        term)
      gsdata[idx, "ymin"] <- i
      gsdata[idx, "ymax"] <- i + 1
      i <- i + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) +
      geom_linerange(aes_(
        ymin = ~ymin,
        ymax = ~ymax,
        color = ~Description
      )) +
      theme_classic(base_size) +
      theme(
        legend.position = "none",
        plot.margin = margin(t = 0, b = 0),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.line.y.left = element_blank(),
        panel.border = element_rect(
          size = 0.25,
          fill = NA,
          colour = "black"
        ),
        plot.title = element_blank(),
        axis.title = element_blank()
      ) +
      scale_x_continuous(expand = c(
        0,
        0
      )) +
      scale_y_continuous(expand = c(0, 0))
    if (length(geneSetID) == 1) {
      v <- seq(1, sum(gsdata$position), length.out = 9)
      inv <- findInterval(rev(cumsum(gsdata$position)), v)
      if (min(inv) == 0) {
        inv <- inv + 1
      }
      col <-
        c(
          rev(RColorBrewer::brewer.pal(5, "Blues")),
          RColorBrewer::brewer.pal(
            5,
            "Reds"
          )
        )
      ymin <- min(p2$data$ymin)
      yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
      xmin <- which(!duplicated(inv))
      xmax <-
        xmin + as.numeric(table(inv)[as.character(unique(inv))])
      d <- data.frame(
        ymin = ymin,
        ymax = yy,
        xmin = xmin,
        xmax = xmax,
        col = col[unique(inv)]
      )
      p2 <- p2 + geom_rect(
        aes_(
          xmin = ~xmin,
          xmax = ~xmax,
          ymin = ~ymin,
          ymax = ~ymax,
          fill = ~ I(col)
        ),
        data = d,
        alpha = 0.9,
        inherit.aes = FALSE
      )
    }
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(
      data = df2,
      aes_(
        x = ~x,
        xend = ~x,
        y = ~y,
        yend = 0
      ),
      color = "grey"
    )
    p.pos <-
      p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") +
      theme(plot.margin = margin(
        t = -0.1,
        r = 0.2,
        b = 0.2,
        l = 0.2,
        unit = "cm"
      ))
    if (!is.null(title) && !is.na(title) && title != "") {
      p.res <- p.res + ggtitle(title)
    }
    if (length(color) == length(geneSetID)) {
      p.res <- p.res + scale_color_manual(values = color)
      if (length(color) == 1) {
        p.res <- p.res + theme(legend.position = "none")
        p2 <- p2 + scale_color_manual(values = "black")
      } else {
        p2 <- p2 + scale_color_manual(values = color)
      }
    }
    if (pvalue_table) {
      pd <- x[geneSetID, c("Description", "NES", "qvalues")]
      colnames(pd)[3] <- c("FDR")
      pd <- pd[, -1]
      pd <- round(pd, 2)
      if (length(geneSetID) == 1) {
        rownames(pd) <- ""
        # pd <- as.data.frame(t(pd))
        # pd$V2 <- c("NES", "FDR")
        # pd <- pd[,c("V2", "V1")]
      }
      print(pd)
      tp <- tableGrob2(pd, p.res)
      p.res <-
        p.res + theme(legend.position = "none") + annotation_custom(
          tp,
          xmin = quantile(p.res$data$x, 0.75),
          xmax = quantile(
            p.res$data$x,
            0.8
          ),
          ymin = quantile(
            p.res$data$runningScore,
            0.6
          ),
          ymax = quantile(
            p.res$data$runningScore,
            0.8
          )
        )
    }
    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(
      axis.line.x = element_line(),
      axis.ticks.x = element_line(),
      axis.text.x = element_text()
    )
    if (length(subplots) == 1) {
      # return(plotlist[[1]] + theme(plot.margin = margin(
      #   t = 0.2,
      #   r = 0.2, b = 0.2, l = 0.2, unit = "cm"
      # )))
      return(plotlist[[1]])
    }
    if (length(rel_heights) > length(subplots)) {
      rel_heights <- rel_heights[subplots]
    }
    aplot::plot_list(
      gglist = plotlist,
      ncol = 1,
      heights = rel_heights
    )
  }

#### 小提琴图 ----
catviolin <- function(adata,
                      x,
                      y,
                      fill = NULL,
                      title = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      legend = TRUE,
                      sort = FALSE,
                      flip = FALSE,
                      hline = NULL,
                      group_by = NULL) {
  x <- enquo(x)
  y <- enquo(y)
  fill <- enquo(fill)
  group_by <- enquo(group_by)

  if (sort) {
    adata <- adata %>%
      mutate({{ x }} := fct_reorder({{ x }}, {{ y }}))
  }
  if (as_label(group_by) != "NULL") {
    adata <- adata %>%
      mutate({{ x }} := fct_reorder({{ x }}, {{ group_by }}))
  }
  p <- adata %>%
    ggplot(aes(x = {{ x }}, y = {{ y }}, fill = {{ fill }})) +
    geom_violin()
  if (flip) {
    p <- p + coord_flip() +
      labs(y = xlab, x = ylab, title = title)
  }
  if (!flip) {
    p <- p +
      labs(x = xlab, y = ylab, title = title)
  }
  p <- p +
    theme_cat() +
    scale_fill_manual(values = cc)
  if (!legend) {
    p <- p + guides(fill = "none")
  }
  if (!is.null(hline)) {
    p <- p + geom_hline(
      yintercept = hline,
      color = "white",
      linetype = "dashed",
      size = 0.5
    )
  }
  return(p)
}


#### scWGCNA网咯 ----
catModuleNetworkPlot <-
  function(seurat_obj,
           mods = "all",
           outdir = "ModuleNetworks",
           plot_size = c(6, 6),
           wgcna_name = NULL,
           label_center = FALSE,
           edge.alpha = 0.25,
           vertex.label.cex = 1,
           vertex.size = 6,
           n_hubs = 25,
           n_center = 10,
           text = "",
           n_conns = 500,
           ...) {
    if (is.null(wgcna_name)) {
      wgcna_name <- seurat_obj@misc$active_wgcna
    }
    MEs <- GetMEs(seurat_obj, wgcna_name)
    modules <- GetModules(seurat_obj, wgcna_name)
    if (mods == "all") {
      mods <- levels(modules$module)
      mods <- mods[mods != "grey"]
    }
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
    cat(paste0("Writing output files to ", outdir))
    TOM <- GetTOM(seurat_obj, wgcna_name)
    n_hubs <- n_hubs
    hub_list <- lapply(mods, function(cur_mod) {
      cur <- subset(modules, module == cur_mod)
      cur <- cur[, c("gene_name", paste0("kME_", cur_mod))] %>%
        top_n(n_hubs)
      colnames(cur)[2] <- "var"
      cur %>%
        arrange(desc(var)) %>%
        .$gene_name
    })
    names(hub_list) <- mods
    for (cur_mod in mods) {
      print(cur_mod)
      cur_color <- modules %>%
        subset(module == cur_mod) %>%
        .$color %>%
        unique()
      n_genes <- n_hubs
      n_conns <- n_conns
      cur_kME <- paste0("kME_", cur_mod)
      cur_genes <- hub_list[[cur_mod]]
      matchind <- match(cur_genes, colnames(TOM))
      reducedTOM <- TOM[matchind, matchind]
      orderind <- order(reducedTOM, decreasing = TRUE)
      connections2keep <- orderind[1:n_conns]
      reducedTOM <- matrix(0, nrow(reducedTOM), ncol(reducedTOM))
      reducedTOM[connections2keep] <- 1
      if (label_center) {
        cur_genes[(n_center + 1):n_hubs] <- ""
      }
      gA <-
        graph.adjacency(
          as.matrix(reducedTOM[1:n_center, 1:n_center]),
          mode = "undirected",
          weighted = TRUE,
          diag = FALSE
        )
      gB <-
        graph.adjacency(
          as.matrix(reducedTOM[
            (n_center + 1):n_genes,
            (n_center + 1):n_genes
          ]),
          mode = "undirected",
          weighted = TRUE,
          diag = FALSE
        )
      layoutCircle <- rbind(layout.circle(gA) / 2, layout.circle(gB))
      g1 <-
        graph.adjacency(
          as.matrix(reducedTOM),
          mode = "undirected",
          weighted = TRUE,
          diag = FALSE
        )
      pdf(
        paste0(outdir, "/", cur_mod, ".pdf"),
        width = plot_size[1],
        height = plot_size[2],
        useDingbats = FALSE
      )
      print(cur_genes)
      print(class(cur_genes))
      cur_genes[!(cur_genes %in% text)] <- ""
      print(cur_genes)
      print(class(g1))
      plot(
        g1,
        edge.color = adjustcolor(cur_color, alpha.f = 0.25),
        edge.alpha = edge.alpha,
        vertex.color = cur_color,
        vertex.label = as.character(cur_genes),
        vertex.label.dist = 1.1,
        vertex.label.degree = -pi / 4,
        vertex.label.color = "black",
        vertex.label.family = "Helvetica",
        vertex.label.font = 3,
        vertex.label.cex = vertex.label.cex,
        vertex.frame.color = "black",
        layout = jitter(layoutCircle),
        vertex.size = vertex.size,
        main = paste(cur_mod)
      )
      dev.off()
    }
  }

#### compass ----
plot_differential_scores <- function(adata, title, text) {
  adata <- adata %>%
    mutate(show_text = if_else(reaction_name %in% text, reaction_name, "")) %>%
    mutate(change = if_else(
      abs(cohens_d) > 0 &
        adjusted_pval <= 0.1,
      if_else(cohens_d > 0, "Up", "Down"),
      "Stable"
    ))
  max_x <- ceiling(max(abs(adata$cohens_d)))

  adata %>% ggplot(aes(
    x = cohens_d,
    y = -log10(adjusted_pval),
    color = change
  )) +
    geom_point() +
    labs(x = "Cohen's d", y = "-log10 (Wilcoxon-adjusted p)", title = title) +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_cat() +
    theme(
      aspect.ratio = 1,
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.margin = margin(t = -10)
    ) +
    scale_color_manual(
      values = c("#00a6e1", "lightgrey", "#ee6470"),
      labels = c(
        str_c(
          "Down (",
          adata %>% count(change) %>% filter(change == "Down") %>% pull(n),
          ")"
        ),
        "Stable",
        str_c(
          "Up (",
          adata %>% count(change) %>% filter(change == "Up") %>% pull(n),
          ")"
        )
      )
    ) +
    ggrepel::geom_text_repel(
      aes(label = show_text),
      max.overlaps = 1000,
      color = "black",
      size = 3
    )
}

plot_all_scores <- function(adata) {
  d <- adata %>%
    filter(adjusted_pval < 0.1) %>%
    group_by(subsystem) %>%
    summarise(d = abs(median(cohens_d))) %>%
    arrange(d) %>%
    pull(subsystem)

  adata %>%
    mutate(alpha = if_else(adjusted_pval <= 0.1, 1, 0.25)) %>%
    mutate(change = if_else(cohens_d > 0, "Up", "Down")) %>%
    mutate(subsystem = factor(subsystem, levels = d)) %>%
    ggplot(aes(
      x = cohens_d,
      y = subsystem,
      alpha = alpha,
      color = change
    )) +
    geom_point() +
    labs(x = "Cohen's d") +
    theme(axis.title.y = element_blank()) +
    scale_color_manual(values = c("#00a6e1", "#ee6470")) +
    theme_cat() +
    guides(color = "none", alpha = "none")
}

#### Seurat ----
catVlnplot <-
  function(seurat_obj,
           features,
           group.by = NULL,
           split.by = NULL,
           sort = FALSE,
           ...) {
    p <- ggrastr::rasterise(
      VlnPlot(
        object = seurat_obj,
        features = features,
        sort = sort,
        pt.size = 0.5,
        group.by = group.by,
        split.by = split.by
      ),
      dpi = 600
    ) +
      NoLegend() + theme_cat() +
      theme(
        axis.title.x = element_blank(),
        plot.title = element_blank(),
        aspect.ratio = 0.5
      )
    p
    return(p)
  }
