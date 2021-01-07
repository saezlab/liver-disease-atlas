# Copyright (c) [2021] [Christian H. Holland]
# christian.holland@bioquant.uni-heidelberg.de

theme_set(cowplot::theme_cowplot())

#' Function that defines my custom theme
#'
#' @param grid character indicating whether grid lines should be plotted. Must
#' be either "x", "y", "no" or NULL.
my_theme <- function(grid = NULL, fsize = 11) {
  if (is.null(grid)) {
    p <- cowplot::background_grid(
      major = "xy", minor = "none",
      size.major = 0.4
    )
  } else if (grid == "y") {
    p <- cowplot::background_grid(
      major = "y", minor = "none",
      size.major = 0.4
    )
  } else if (grid == "x") {
    p <- cowplot::background_grid(
      major = "x", minor = "none",
      size.major = 0.4
    )
  } else if (grid == "no") {
    p <- theme()
  }
  p <- p +
    theme(
      title = element_text(size = fsize + 1),
      axis.text = element_text(size = fsize),
      legend.text = element_text(size = fsize),
      legend.title = element_text(size = fsize),
      plot.subtitle = element_text(size = fsize),
      strip.background = element_rect(colour = "white", fill = "white"),
      strip.text = element_text(size = fsize)
    )
}


#' This functions plots a PCA plots
plot_pca <- function(pca_result, x = "PC1", y = "PC2", feature = NULL) {
  if (!is.null(feature)) {
    if (!(feature %in% colnames(pca_result$coords))) {
      stop(feature, " is not a valid feature")
    }

    p <- pca_result$coords %>%
      ggplot(aes(
        x = .data[[as.character(x)]],
        y = .data[[as.character(y)]],
        color = .data[[as.character(feature)]]
      ))
  } else if (is.null(feature)) {
    p <- pca_result$coords %>%
      ggplot(aes(
        x = .data[[as.character(x)]],
        y = .data[[as.character(y)]]
      ))
  }
  p +
    geom_point(size = 2) +
    labs(
      x = paste0(x, " (", pca_result$var[parse_number(x)], "%)"),
      y = paste0(y, " (", pca_result$var[parse_number(y)], "%)"),
      color = feature
    )
}

#' This functions plots the library size of RNA-seq samples as barplot.
#'
#' @param count_matrix count matrix with genes in rows and samples in columns.
plot_libsize <- function(count_matrix) {
  df <- colSums(count_matrix) %>%
    enframe("sample", "libsize")

  ggplot(df, aes(x = libsize, y = fct_reorder(sample, libsize))) +
    geom_col() +
    labs(x = "Library size", y = "Sample")
}


#' This function plots volcano plot(s)
#'
#' A volcano plot it plotted for each contrast. In case that multiple contrasts
#'   are provided each contrast will be an individual facet.
#'
#' @param df A data frame that stores the output of a differential gene
#'   expression analysis (e.g. via \code{\link{limma}}). This data frame must
#'   contain the columns \code{logFC, contrast, regulation} and one of
#'   \code{pval, p, p-value, p.value}. The column \code{regulation} should be a
#'   factor containing three levels reflecting up-, down and non significant
#'   regulation. The column \code{contrast} should contain an identifier for
#'   this contrast. It is possible to have multiple contrast within the data
#'   frame. All other columns will be ignored.
#' @param ... further parameters of \code{facet_wrap()} function, e.g.
#'   \code{scales="free"}
#'
#' @return ggplot object of volcano plots
plot_volcano <- function(df, ...) {
  df %>%
    rename(p = any_of(c("pval", "p", "p-value", "p.value"))) %>%
    ggplot(aes(
      x = logFC, y = -log10(p), color = regulation,
      alpha = regulation
    )) +
    geom_point() +
    lemon::facet_rep_wrap(~contrast, ...) +
    labs(x = "logFC", y = expression(-log["10"] * "(p-value)")) +
    scale_color_manual(values = AachenColorPalette::aachen_color(
      c("red", "blue", "black50")
    ), drop = F) +
    scale_alpha_manual(values = c(0.7, 0.7, 0.2), guide = "none", drop = F)
}


#' This function plots p-value histogram(s)
#'
#' @param A data frame that must contain at least a column with p-values. The
#'   column must be named either: \code{pval, p, p-value, p.value}. If there
#'   exist different groups (e.g. in case of GSEA on different contrast each
#'   contrast is a group) their identifier should be stored in an additional
#'   column with variable column name. All other columns will be ignored
#' @param facet_var Quoted string of the column name that stores the group
#'   information. A p-value histogram will be plotted for each group. In case
#'   of multiple groups also a vector of quoted strings can be provided.
#' @param ... further parameters of \code{facet_wrap()} function, e.g.
#'   \code{scales="free"}
#'
#' @return ggplot object of p-value histograms
plot_phist <- function(df, facet_var = NULL, ...) {
  p_synonyms <- c("pval", "p", "p-value", "p.value", "p_val", "P.Value")

  p <- df %>%
    rename(p = any_of(p_synonyms)) %>%
    ggplot(aes(x = p)) +
    geom_histogram(color = "white", bins = 21, boundary = 0) +
    labs(x = "p-value", y = "Count") +
    scale_x_continuous(limits = c(0, 1.0001), labels = label_percent())

  if (!is.null(facet_var)) {
    if (length(facet_var) > 1) {
      facet_var <- str_c(facet_var, collapse = "+")
    }
    p <- p +
      lemon::facet_rep_wrap(as.formula(str_c("~", facet_var)), ...)
  }
  return(p)
}

#' Venn diagram
#'
#' This function plots a venn diagram summarising the overlap of differential
#'   expressed genes separately for up and down regulated genes across various
#'   classes.
#'
#' @param tables list of tables. Each individual table must contain the columns
#'   \code{gene}, \code{class} and \code{regulation}. Each set in the final plot
#'   is a combination of the column \code{class} and \code{regulation}. All
#'   other columns will be ignored.
#' @param whether logical indicating whether the venn diagrams should be merged
#'   as a collage
#'
#' @return Venn diagram
plot_venn_diagram <- function(tables = list(), collage = TRUE) {
  if (length(tables) == 2) {
    c1 <- tables[[1]] %>%
      distinct(class) %>%
      pull() %>%
      as.character()
    c2 <- tables[[2]] %>%
      distinct(class) %>%
      pull() %>%
      as.character()
    t1 <- tables[[1]] %>% count(regulation)
    t2 <- tables[[2]] %>% count(regulation)

    plots <- c("up", "down") %>%
      map(function(r) {
        # set sizes of regulated genes
        a1 <- t1 %>%
          filter(regulation == r) %>%
          pull(n)
        a2 <- t2 %>%
          filter(regulation == r) %>%
          pull(n)
        ca <- intersect(
          tables[[1]] %>% filter(regulation == r) %>% pull(gene),
          tables[[2]] %>% filter(regulation == r) %>% pull(gene)
        ) %>%
          length()


        grid.newpage()
        draw.pairwise.venn(
          area1 = a1, area2 = a2, cross.area = ca,
          category = c(c1, c2),
          lty = "blank",
          cex = 1,
          fontfamily = rep("sans", 3),
          fill = AachenColorPalette::aachen_color(c("purple", "petrol")),
          cat.col = AachenColorPalette::aachen_color(c("purple", "petrol")),
          cat.cex = 1.1,
          cat.fontfamily = rep("sans", 2)
        ) %>%
          as_ggplot() %>%
          grid.arrange(top = textGrob(r, gp = gpar(
            fontsize = 15,
            fontface = "bold"
          )))
      })
  } else if (length(tables) == 3) {
    c1 <- tables[[1]] %>%
      distinct(class) %>%
      pull() %>%
      as.character()
    c2 <- tables[[2]] %>%
      distinct(class) %>%
      pull() %>%
      as.character()
    c3 <- tables[[3]] %>%
      distinct(class) %>%
      pull() %>%
      as.character()
    t1 <- tables[[1]] %>% count(regulation)
    t2 <- tables[[2]] %>% count(regulation)
    t3 <- tables[[3]] %>% count(regulation)

    plots <- c("up", "down") %>%
      map(function(r) {
        # set sizes of regulated genes
        a1 <- t1 %>%
          filter(regulation == r) %>%
          pull(n)
        a2 <- t2 %>%
          filter(regulation == r) %>%
          pull(n)
        a3 <- t3 %>%
          filter(regulation == r) %>%
          pull(n)

        a12 <- purrr::reduce(
          list(
            tables[[1]] %>% filter(regulation == r) %>% pull(gene),
            tables[[2]] %>% filter(regulation == r) %>% pull(gene)
          ),
          intersect
        ) %>%
          length()

        a23 <- purrr::reduce(
          list(
            tables[[2]] %>% filter(regulation == r) %>% pull(gene),
            tables[[3]] %>% filter(regulation == r) %>% pull(gene)
          ),
          intersect
        ) %>%
          length()

        a13 <- purrr::reduce(
          list(
            tables[[1]] %>% filter(regulation == r) %>% pull(gene),
            tables[[3]] %>% filter(regulation == r) %>% pull(gene)
          ),
          intersect
        ) %>%
          length()

        a123 <- purrr::reduce(
          list(
            tables[[1]] %>% filter(regulation == r) %>% pull(gene),
            tables[[2]] %>% filter(regulation == r) %>% pull(gene),
            tables[[3]] %>% filter(regulation == r) %>% pull(gene)
          ),
          intersect
        ) %>%
          length()

        grid.newpage()
        p <- draw.triple.venn(
          area1 = a1, area2 = a2, area3 = a3,
          n12 = a12, n23 = a23, n13 = a13,
          n123 = a123,
          category = c(c1, c2, c3),
          lty = "blank",
          cex = 1,
          fontfamily = rep("sans", 7),
          fill = AachenColorPalette::aachen_color(c(
            "purple", "petrol",
            "red"
          )),
          cat.col = AachenColorPalette::aachen_color(c(
            "purple", "petrol",
            "red"
          )),
          cat.cex = 1.1,
          cat.fontfamily = rep("sans", 3)
        ) %>%
          as_ggplot() %>%
          grid.arrange(top = textGrob(r, gp = gpar(
            fontsize = 15,
            fontface = "bold"
          )))
      })
  }

  if (collage) {
    pp <- plot_grid(plotlist = plots)
    return(pp)
  } else {
    return(plots)
  }
}

#' Upset plot
#'
#' This function plots an upset plot summarizing the overlap of differential
#'   expressed genes across various classes.
#'
#' @parram tables list of tables. Each individual table must contain the columns
#'   \code{gene}, \code{class} and \code{regulation}. Each set in the final plot
#'   is a combination of the column \code{class} and \code{regulation}. All
#'   other columns will be ignored.
#' @params A list of named options to pass to
#' \code{\link[=upset]{UpSetR::upset()}} such as \code{nsets} or
#'   \code{nintersects}. These options should not \code{mat}, \code{text.scale}
#'   or \code{point.size}.
#'
#' @return An upset plot
plot_upset <- function(tables = list(),
                       options = list(), ...) {
  df <- do.call(bind_rows, tables)

  mat <- df %>%
    filter(regulation != "ns") %>%
    unite(key, class, regulation) %>%
    select(gene, key) %>%
    mutate(val = 1) %>%
    spread(key, val, fill = 0) %>%
    data.frame(row.names = 1, stringsAsFactors = F, check.names = F)

  plot <- do.call(
    upset,
    c(
      list(
        data = mat,
        mainbar.y.label = "Common genes",
        sets.x.label = "Total number of genes",
        text.scale = 2, point.size = 3
      ),
      options
    )
  )

  return(plot)
}



#' Plot of temporal expression profiles from stem analysis
#'
#' This function plots all significant temporal expression profile cluster as
#'   facets. Each facet contains the expression of all gene members and a
#'   summarized trajectory across all gene members.
#'
#' @param df Output from \code{run_stem()} function. In case that the output
#'   contains multiple stem runs the output must be filtered beforehand for a
#'   single stem run
#' @param p_cutoff numeric value setting the cutoff of significant temporal
#'   expression profiles
#' @param model_profile logical indicating whether the profile models from stem
#'   should be plotted along the individual and summary trajectories.
#'
#' @return ggplot of temporal expression profiles per cluster
plot_stem_profiles <- function(df, p_cutoff = 0.05, model_profile = F,
                               features = "genes", ...) {
  # ensure that only one key is present in input
  if (n_distinct(df$key) != 1) {
    stop(
      "Only one 'key' must be present. Currenty the keys are: ",
      str_c(unique(df$key), collapse = ", ")
    )
  }

  # # order profiles by significance
  # df2 = df %>%
  #   filter(p <= p_cutoff) %>%
  #   arrange(p) %>%
  #   mutate(profile = fct_inorder(as.character(profile))) %>%
  #   group_by(profile) %>%
  #   mutate(rank = cur_group_id()) %>%
  #   ungroup() %>%
  #   mutate(profile = str_c("Cluster ", rank))

  # get mean profile coords
  profile_summary <- df %>%
    group_by(key, profile, p, time) %>%
    dplyr::summarize(
      m = median(value),
      s = sd(value),
      lower = m - s,
      upper = m + s
    ) %>%
    ungroup()


  # # extract meta data of profiles
  # profile_anno = df %>%
  #   group_by(key, profile, p) %>%
  #   mutate(y = 1.2*abs(max(value))) %>%
  #   ungroup() %>%
  #   mutate(max_time = max(time)) %>%
  #   distinct(key, profile, p, size, y, max_time) %>%
  #   # mutate(label = str_c(size, " ", features, "\n(p = ", p, ")"))
  #   mutate(label = str_c(size, " ", features))

  p <- df %>%
    ggplot(aes(x = time, y = value, group = gene)) +
    geom_hline(yintercept = 0) +
    geom_line(alpha = 0.05) +
    geom_line(
      data = profile_summary, aes(x = time, y = m), inherit.aes = F,
      color = AachenColorPalette::aachen_color("green"), size = 0.5
    ) +
    geom_errorbar(
      data = profile_summary, aes(
        x = time, y = m, ymin = lower,
        ymax = upper
      ),
      inherit.aes = F,
      color = AachenColorPalette::aachen_color("green"),
      size = 0.55, width = 0.25
    ) +
    lemon::facet_rep_wrap(~profile, scales = "free", ...) +
    # my_theme(grid = "no") +
    scale_x_continuous(breaks = unique(df$time)) +
    labs(x = "Time", y = "logFC") +
    theme(strip.text = element_text(angle = 0, hjust = 0))

  if (model_profile) {
    # get coords of stem model profiles
    profile_stem_summary <- df %>%
      distinct(key, profile, y_coords) %>%
      separate_rows(y_coords, sep = ",") %>%
      mutate(y_coords = as.numeric(y_coords)) %>%
      group_by(key, profile) %>%
      mutate(time = unique(df$time)) %>%
      ungroup()

    p <- p +
      geom_line(
        data = profile_stem_summary, aes(x = time, y = y_coords),
        inherit.aes = F,
        color = AachenColorPalette::aachen_color("blue"), size = 1
      )
  }
  return(p)
}


#' Plot a Gene Set Enrichment Analysis (GSEA) plot for a given gene signature
#' and gene set
#'
#' @param signature dataframe of a gene signature. The dataframe must contain at
#'  least the gene ids in a column named \code{gene} and a column (with custom
#'  column name) containing gene level statistics (e.g. logFC or t-values).
#'  Additional columns will be ignored.
#' @param geneset dataframe of a gene set. The dataframe must contain at
#'  least the gene ids in a column named \code{gene}. Additional columns will be
#'  ignored. Note that both dataframes \code{signature} and \code{geneset}
#'  must contain the same gene id type (e.g. gene symbols)
#' @param gene_level_stat unquoted string indicating the name of the gene
#' level statistic column (that is stored in the dataframe \code{signature}).
#'
#' @return A GSEA-plot (ggplot object)
make_gsea_plot <- function(signature, geneset, gene_level_stat, ...) {

  # scaling and ranking (decreasing) gene level statistic (e.g. logFC)
  stats_df <- signature %>%
    rename(stat := {{ gene_level_stat }}) %>%
    # rename(stat = logFC) %>% # uncommend line to test the function
    arrange(-stat) %>%
    transmute(gene,
      stat = stat / max(abs(stat)),
      rank = row_number()
    )

  # convert datafrane to named list of gene ids and scaled gene level statistic
  stats_list <- stats_df %>%
    select(gene, stat) %>%
    deframe()

  # extract the ranks/indices of the geneset member within the ranked signature
  indices_of_geneset_member <- stats_df %>%
    inner_join(geneset, by = "gene") %>%
    arrange(rank) %>%
    pull(rank)

  # run gsea
  gseaRes <- calcGseaStat(stats_list,
    selectedStats = indices_of_geneset_member,
    returnAllExtremes = TRUE
  )

  # extrat results from gsea analysis
  bottoms <- gseaRes %>% pluck("bottoms")
  tops <- gseaRes %>% pluck("tops")

  # highest rank = number of genes in signature
  max_rank <- stats_df %>%
    pull(rank) %>%
    max()

  # construct coordinates of enrichment plot
  xs <- rep(indices_of_geneset_member, each = 2)
  ys <- rbind(bottoms, tops) %>% as.vector()

  # build plotting table
  enrichment_df <- tibble(
    rank = c(0, xs, max_rank + 1),
    running_sum = c(0, ys, 0),
    max_top = max(tops),
    min_bottom = min(bottoms)
  )

  # position of the geneset member in the ranked signature
  gene_positions <- tibble(rank = indices_of_geneset_member)

  # extraction of highest absolute enrichment score and corresponding rank
  top_es <- enrichment_df %>%
    top_n(1, abs(running_sum))

  # plotting of enrichment curve
  p1 <- ggplot(data = enrichment_df, aes(x = rank, y = running_sum)) +
    geom_hline(aes(yintercept = max_top),
      colour = "#A11035", linetype = "dashed"
    ) +
    geom_hline(aes(yintercept = min_bottom),
      colour = "#A11035", linetype = "dashed"
    ) +
    geom_hline(yintercept = 0, colour = "black") +
    geom_segment(
      data = top_es,
      mapping = aes(x = rank, y = 0, xend = rank, yend = running_sum),
      linetype = "dashed"
    ) +
    geom_path(size = 1, color = "#0098A1") +
    lims(x = c(0, max_rank + 1)) +
    labs(y = "Enrichment Score") +
    background_grid(major = "y", minor = "none", size.major = 0.4) +
    theme(
      axis.line.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = unit(c(1, 1, -0.25, 1), "cm"),
      legend.position = "none",
      title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    ) +
    annotate("text",
      x = 1, y = max(tops) + 0.05,
      label = str_c("ES:", round(gseaRes$res, 2), sep = " "),
      size = 4.5, hjust = 0
    )

  # plotting of gene positions among the signature
  p2 <- ggplot(stats_df, aes(x = rank, y = 1)) +
    geom_tile(aes(color = rank)) +
    geom_segment(gene_positions,
      mapping = aes(x = rank, y = 1.51, xend = rank, yend = 3.51),
      size = 0.5, color = "black", alpha = 0.5
    ) +
    scale_color_gradientn(
      colours = c(
        "#A11035", "#F6A800", "#FFED00",
        "#57AB27", "#00549F", "#612158"
      ),
      breaks = c(1, max_rank / 2, max_rank),
      labels = c("High", "", "Low")
    ) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(),
      plot.margin = unit(c(-0.25, 1, 1, 1), "cm"),
      axis.line.x = element_blank(),
      title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.position = "bottom",
      legend.justification = "center"
    ) +
    labs(x = "Rank", y = NULL, color = "Gene-level Statistic") +
    guides(color = guide_colorbar(barwidth = 15, ticks = F, title.vjust = 0.8, title.position = "top"))

  # combine both plots
  p <- plot_grid(p1, p2, ncol = 1, align = "v", axis = "l", rel_heights = c(2, 1))
  return(p)
}

#' This funtion plots top differentially expressed genes from a differential
#' expression analysis.
plot_top_genes <- function(data, class, fontsize, ...) {
  max_logfc <- data %>%
    slice_max(order_by = abs(logFC), n = 1) %>%
    pull(logFC) %>%
    abs()

  data <- data %>%
    mutate(stars = gtools::stars.pval(fdr))

  up_plot <- data %>%
    filter(logFC >= 0) %>%
    ggplot(aes(x = logFC, y = fct_reorder(gene, logFC), fill = as_factor(sign(logFC)))) +
    geom_col() +
    labs(x = "Up-regulated logFC", y = NULL) +
    scale_x_continuous(
      position = "top", limits = c(0, max_logfc),
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))
    ) +
    theme(legend.position = "none") +
    scale_fill_manual(values = AachenColorPalette::aachen_color("red")) +
    my_theme(grid = "no", fsize = fz) +
    geom_text(aes(label = stars),
      vjust = +0.75, hjust = 1.5,
      color = "white", size = fontsize / (14 / 5)
    ) +
    theme(plot.margin = margin(0, -0.15, 0, 0, "cm"))

  down_plot <- data %>%
    filter(logFC < 0) %>%
    ggplot(aes(x = -logFC, y = fct_reorder(gene, logFC), fill = as_factor(sign(logFC)))) +
    geom_col() +
    labs(x = "Down-regulated logFC", y = NULL) +
    scale_x_reverse(
      limits = c(max_logfc, 0),
      breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))
    ) +
    scale_y_discrete(position = "right") +
    theme(legend.position = "none") +
    scale_fill_manual(values = AachenColorPalette::aachen_color("blue")) +
    my_theme(grid = "no", fsize = fz) +
    geom_text(aes(label = stars),
      vjust = +0.75, hjust = -0.5,
      color = "white", size = fontsize / (14 / 5)
    ) +
    theme(plot.margin = margin(0, 0, 0, -0.15, "cm"))

  up_plot + down_plot
}

#' This functions plots a density distribution of GO terms of clusters.
#'
#' @param tibble containing the columns `rank` and `description` and
#' `regulation`.
#'
#' @return density plot.
plot_go_rank_density <- function(df) {
  df %>%
    ggplot(aes(x = rank, colour = description)) +
    stat_density(geom = "line", position = "identity") +
    geom_rug() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "top"
    ) +
    labs(
      x = str_wrap("GO term ranking based on p-value"), y = "Density",
      color = NULL
    ) +
    lemon::facet_rep_wrap(~regulation, scales = "free", drop = T) +
    scale_x_continuous(breaks = function(x) {
      as.numeric(
        gsub("^0", 1, unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))))
      )
    })
}

#' This functions plots a wordcloud
#'
#' @param df tibble with two columns `word` and `n`. Other columns will be
#' ignored.
#' @fontsize integer indicating desired fontsize
#'
#' @return wordcloud
plot_wordcloud <- function(df, fontsize = 9) {
  df %>%
    ggplot(aes(label = word, size = n)) +
    geom_text_wordcloud() +
    scale_size_area(max_size = fontsize / (14 / 5)) +
    theme(axis.line = element_blank())
}
