# Copyright (c) [2021] [Christian H. Holland]
# christian.holland@bioquant.uni-heidelberg.de

#' This functions performs PCA analysis either on all genes or only on the most
#' variable genes.
do_pca <- function(df, meta = NULL, top_n_var_genes = NULL) {
  if (!is.null(top_n_var_genes)) {
    top_var_genes <- base::apply(df, 1, var) %>%
      sort() %>%
      tail(top_n_var_genes) %>%
      names()

    df <- df[top_var_genes, , drop = F]
  }
  pca_obj <- prcomp(t(df), center = T, scale. = T)

  coords <- pca_obj$x %>%
    data.frame(check.names = F, stringsAsFactors = F) %>%
    rownames_to_column("sample") %>%
    as_tibble()

  if (!is.null(meta)) {
    coords <- left_join(coords, meta, by = "sample")
  }

  var <- round(summary(pca_obj)$importance[2, ] * 100, 2)

  res <- list()
  res$coords <- coords
  res$var <- var

  return(res)
}

#' This function performs differential gene expression analysis via the limma
#'   workflow.
#'
#' @param expr expression matrix with genes in rows and samples in columns.
#' @param design design matrix. For more information see
#'   \code{\link[limma:lmFit]{lmFit}}.
#' @param contrasts contrasts. For more information see
#'   \code{\link[limma:contrasts.fit]{contrasts.fit}}. If NULL no contrasts are
#'   computed
#'
#' @return limma result in table format. Contains the columns \code{gene,
#'   (contrast), logFC, statistic, pval, fdr}.
run_limma <- function(expr, design, contrasts = NULL, ...) {
  if (!is.null(contrasts)) {
    limma_result <- lmFit(expr, design) %>%
      contrasts.fit(contrasts) %>%
      eBayes() %>%
      tidy() %>%
      select(gene,
        contrast = term, logFC = estimate, statistic = statistic,
        pval = p.value
      ) %>%
      group_by(contrast) %>%
      mutate(fdr = p.adjust(pval, method = "BH")) %>%
      ungroup()
  } else if (is.null(contrasts)) {
    limma_result <- lmFit(expr, design) %>%
      eBayes() %>%
      topTableF(n = Inf) %>%
      rownames_to_column("gene") %>%
      as_tibble()
  }

  return(limma_result)
}

#' Classification/Assignment of differential expressed genes
#'
#' This function classifies genes as deregulated based on effect size and
#'   adjusted p-value. A gene is considered as differential expressed if the
#'   adjusted p-value is below AND if the absolute effect size is above user
#'   specified cutoffs.
#'
#' @param df A table that must contain at least the columns \code{gene},
#'   \code{fdr} and \code{logFC}. The table must not include a column named
#'   \code{regulation}.
#' @param fdr_cutoff numeric value that denotes the adjusted p-value cutoff.
#' @param effect_size_cutoff numeric value that denotes the effect size cutoff.
#'
#' @return The input table with an additional colum names \code{regulation}. A
#'   gene can be upregulated (up), downregulated (down) or not significantly
#'   regulated (ns).
assign_deg <- function(df, fdr_cutoff = 0.05, effect_size_cutoff = 1,
                       fdr_id = fdr, effect_size_id = logFC) {
  degs <- df %>%
    mutate(regulation = case_when(
      {{ effect_size_id }} >= effect_size_cutoff & {{ fdr_id }} <= fdr_cutoff ~
      "up",
      {{ effect_size_id }} <= -effect_size_cutoff & {{ fdr_id }} <= fdr_cutoff ~
      "down",
      TRUE ~ "ns"
    )) %>%
    mutate(regulation = factor(regulation, levels = c("up", "down", "ns")))

  return(degs)
}

#' Tidy a matrix
#'
#' This utility function takes a matrix and converts it to a tidy format and
#' adds if available observations' meta data.
#'
#' @param mat A matrix with observations/features in rows and variables in
#' columns
#' @param feature Class name of observations/features, e.g.
#' transcription_factors
#' @param key Class name of variables, e.g. samples
#' @param value Class name of matrix values, e.g. activities
#' @param meta Data frame with meta data of the observations. To map the meta
#' data to the tidied table the observation/feature column name must be
#' identical.
#'
#' @return Tidy table.
#'
#' @export
tdy <- function(mat, feature, key, value, meta = NULL) {
  res <- mat %>%
    data.frame(check.names = F, stringsAsFactors = F) %>%
    rownames_to_column(feature) %>%
    as_tibble() %>%
    gather({{ key }}, {{ value }}, -{{ feature }})

  if (!is.null(meta)) {
    res <- res %>%
      left_join(meta, by = key)
  }

  return(res)
}

#' Untidy a tibble
#'
#' This utility function takes a tidy tibble and converts it to a matrix.
#'
#' @param tbl A tidy tibble
#' @param feature Class name of observations/features present in tidy tibble
#' @param key Class name of key present in tidy tibble
#' @param value Class name of values in tidy tibble
#'
#' @return Matrix with observation in rows and variables in columns.
#'
#' @export
untdy <- function(tbl, feature, key, value, fill = NA) {
  tbl %>%
    select({{ feature }}, {{ key }}, {{ value }}) %>%
    spread({{ key }}, {{ value }}, fill = fill) %>%
    data.frame(row.names = 1, check.names = F, stringsAsFactors = F)
}


#' Make gene sets for GSEA
#'
#' This function convert gene sets in a table format to the format required by
#' the \code{\link[=viper]{fgsea::fgsea}} function.
#'
#' @param genesets A dataframe of gene set that must contain the
#' columns \code{geneset} and \code{gene}.
#'
#' @return Gene sets in the \code{\link[=fgsea]{fgsea::fgsea()}} format.
#'
#' @export
make_gsea_genesets <- function(genesets) {
  genesets %>%
    group_by(geneset) %>%
    summarise(tmp = list(gene)) %>%
    deframe()
}


#' GSEA wrapper
#'
#' This function is a convenient wrapper for the
#' \code{\link[=fgsea]{fgsea::fgsea()}} function to perform GSEA
#'
#' @param sig_df An expression matrix with genes (HGNC symbol) in rows and
#' signatures in columns. Signatures can be either single samples or contrats.
#' Potential NA's won't be considered.
#' @param genesets A data frame of gene sets that must contain at least the
#' columns \code{geneset} and \code{gene}.
#' @param nperm Integer indicating the number of permutations to generate the
#' null distribution.
#' @param options A list of named options to pass to
#' \code{\link[=fgsea]{fgsea::fgsea()}} such as  or \code{minSize}. These
#' options should not include \code{pathways}, \code{stats} or \code{nperm}.
#' @param tidy Logical, whether computed normalized enrichment scores should be
#' returned in a tidy format.
#'
#' @return A matrix of normalized enrichment score for each gene sets across all
#'  signatures. If \code{tidy} is TRUE the results are retured in a tidy format,
#'  including also p-values and leading edge genes.
run_gsea <- function(sig_df, genesets, nperm = 1000, options = list(),
                     tidy = F) {
  geneset_list <- make_gsea_genesets(genesets)
  gsea_res <- apply(sig_df, 2, function(col) {
    do.call(
      fgsea,
      c(
        list(
          pathways = geneset_list,
          stats = col[!is.na(col)], nperm = nperm
        ),
        options
      )
    )
  }) %>%
    enframe("signature", "value") %>%
    unnest(value) %>%
    rename(geneset = pathway)

  if (tidy) {
    metadata <- genesets %>%
      select(-c(gene)) %>%
      distinct()

    tidy_gsea_res <- left_join(gsea_res, metadata, by = "geneset")
    return(gsea_res)
  } else {
    return(untdy(gsea_res, "geneset", "signature", "NES"))
  }
}


#' This function performs an over-representation analysis using Fishers exact
#' test
#'
#' @param sig a datframe containing a signature you would like to analyse with
#'  various gene sets. This dataframe must contain at least the column
#'  \code{gene} containing the genes of the signature to analyse. All other
#'  columns will be ignored.
#' @param sets A data frame holding all gene sets. It must contain at least the
#'   columns \code{geneset}, \code{gene} and \code{group}, containing gene set
#'   identifier, gene set members and a variable separating gene sets into
#'   groups, respectively (e.g. GO and KEGG). Further meta information of the
#'   gene sets might be provided in further columns.
#' @param min_size integer indication the minimum gene set size. Gene sets with
#'   less members are discarded from the analysis.
#' @param options A list of named options to pass to
#'   \code{\link[=fisher.test]{base::fisher.test()}} such as  or
#'   \code{alternative}.
#' @param background_n Integer indicating the background size of the signature.
#'   If not specified the number of background genes it determined by the total
#'   number of unique genes in the gene sets.
#'
#' @return Table reporting for each gene set the contingency table, estimate,
#'   side of the test, p-value and fdr.
run_ora <- function(sig, sets, min_size = 10, options = list(),
                    background_n = NULL) {

  #'                       signature (sig)
  #'                     yes   |   no    ||   total
  #'                   ------------------------------
  #'           yes   |     2   |   117   || 119
  #'                 | ------------------------------
  #'     set   no    |     72  |   4389  || 4461
  #'                 |===============================
  #'           total |     74  |   4506  || 4580


  # exclude gene sets that are too small
  sets <- sets %>%
    add_count(geneset, group) %>%
    filter(n >= min_size) %>%
    select(-n)

  # determine background_n if not specified by the user
  if (is.null(background_n)) {
    background_n <- length(unique(sets$gene))
  }

  res <- sets %>%
    nest(set = c(gene)) %>%
    # construction of contingency table
    mutate(contingency_table = set %>% map(function(set) {
      sig_vec <- sig$gene
      set_vec <- set$gene

      sig_total <- length(sig_vec)
      set_total <- length(set_vec)

      overlap_n <- intersect(sig_vec, set_vec) %>% length()
      sig_specific_n <- setdiff(sig_vec, set_vec) %>% length()
      set_specific_n <- setdiff(set_vec, sig_vec) %>% length()

      no_no_n <- background_n - sig_total - set_specific_n


      t <- matrix(c(overlap_n, sig_specific_n, set_specific_n, no_no_n),
        nrow = 2, dimnames = list(
          Set = c("yes", "no"),
          Signature = c("yes", "no")
        )
      )

      return(t)
    })) %>%
    # run fisher test
    mutate(stat = contingency_table %>% map(function(ct) {
      if (ct[1, 1] == 0) {
        stat <- as_tibble(NULL)
      } else {
        stat <- do.call(fisher.test, c(
          list(x = ct),
          options
        )) %>%
          tidy()
      }

      return(stat)
    })) %>%
    unnest(stat) %>%
    group_by(group) %>%
    mutate(fdr = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    select(-set, -conf.low, -conf.high, -method)

  return(res)
}

#' This function translates genes across different gene id classes from mouse
#'   and human
#'
#' @param x Data frame that must store at least a column containing genes.
#' @param from string indicating the present gene id class. Must match to a
#'   column in the annotation data frame (e.g.
#'   \code{symbol_mgi, ensembl_mgi, entrez_hgnc, ...})
#' @param to string indicating the desired gene id class. Similar to
#'   \code{from} parameter.
#'
#' @param ref string indicating the column name of the genes. Mostly \code{gene}
#' @param na_rm gene id mapping can lead to NAs. This logical argument indicates
#'   whether rows with NA should be removed.
#'
#' @return The input is returned with translated gene ids.
translate_gene_ids <- function(x, from, to, ref = "gene", na_rm = F) {
  anno <- readRDS(here::here("data/annotation/gene_id_annotation.rds"))

  if (!all(c(from, to) %in% colnames(anno))) {
    stop(
      "'from' and 'to' must be one of: ",
      str_c(colnames(anno), collapse = ", ")
    )
  }

  annotation_df <- anno %>%
    select({{ from }}, {{ to }}) %>%
    drop_na()

  mapped_x <- x %>%
    rename({{ from }} := {{ ref }}) %>%
    left_join(annotation_df, by = from) %>%
    select({{ ref }} := {{ to }}, everything(), -{{ from }})

  if (na_rm == T) {
    res <- mapped_x %>%
      drop_na({{ ref }})
  } else {
    res <- mapped_x
  }
  return(res)
}


#' This function computes the similarity of feature sets
#'
#' @param mat A dataframe with feature set ids in columns. The matrix is filled
#'   with feature set members. Row names should be a simplex index (1,2,3,...).
#'   In case of unbalanced feature set size the matrix should be filled up with
#'   NA that will be ignored for the similarity computation.
#' @param measure String selecting which measure should be used to compute
#'   feature set similarity. Either Jaccard index (jaccard) or Overlap
#'   Coefficient (overlap_coef) can be chosen
#' @tidy logical indicating whether this function should return the similarity
#'   matrix in a tidy format
#'
#' @return A (tidy) similarity matrix denoting the pairwise similarity between
#'    all provided feature sets
set_similarity <- function(mat, measure = "jaccard", tidy = T) {

  # number of columns
  n <- ncol(mat)

  # set up similarity matrix
  similarity_mat <- matrix(NA, n, n)
  colnames(similarity_mat) <- colnames(mat)
  rownames(similarity_mat) <- colnames(mat)

  # set the diagonal to 0
  diag(similarity_mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      s1 <- na.omit(mat[, i])
      s2 <- na.omit(mat[, j])

      if (measure == "jaccard") {
        jaccard_index <- length(intersect(s1, s2)) / length(union(s1, s2))
        similarity_mat[i, j] <- jaccard_index
      } else if (measure == "overlap_coef") {
        overlap_coef <- length(intersect(s1, s2)) / min(length(s1), length(s2))
        similarity_mat[i, j] <- overlap_coef
      } else if (measure == "intersect") {
        common <- length(intersect(s1, s2))
        s1_specific <- length(setdiff(s1, s2))
        s2_specific <- length(setdiff(s2, s1))

        similarity_mat[i, j] <- str_c(common, s1_specific, s2_specific,
          sep = "|"
        )
      }
    }
  }

  if (tidy == T & measure != "intersect") {
    res <- similarity_mat %>%
      t() %>%
      data.frame(check.names = F, stringsAsFactors = F) %>%
      tdy("set1", "set2", "similarity") %>%
      drop_na() %>%
      filter(set1 != set2) %>%
      mutate_if(is.character, factor, levels = colnames(mat)) %>%
      mutate(set1 = fct_rev(set1))
  } else if (tidy == T & measure == "intersect") {
    res <- similarity_mat %>%
      t() %>%
      data.frame(check.names = F, stringsAsFactors = F) %>%
      tdy("set2", "set1", "similarity") %>%
      drop_na() %>%
      filter(set1 != set2) %>%
      select(set1, set2, similarity) %>%
      separate(similarity,
        into = c("common", "set1_specific", "set2_specific"),
        convert = T
      ) %>%
      mutate_if(is.character, factor, levels = colnames(mat)) %>%
      mutate(set1 = fct_rev(set1))
  } else {
    res <- similarity_mat
  }

  return(res)
}



#' This functions clusters temporal expression profiles using the STEM software
#'
#' @param path path to a directory that itself must contain the directories
#'   \code{input, setting, output}. The directory \code{input} contains
#'   expression files that sever as input for stem analyses. The directory
#'   \code{setting} contains (muliple) files storing parameters of the stem
#'   analysis. The parameter \code{Data_File} must match an input file in the
#'   directory \code{input}, including path relative to working directory. It
#'   is possible to store multiple setting files in  this directory, stem will
#'   perform an analysis for all of them. The directory \code{output} will be
#'   used to store the results/output of the stem analysis.
#' @param jar_path path to the java executable (stem.jar) relative to the
#'   working directory
#' @param clear_output logical indicating whether the (intermediate) output
#'   files of stem should be removed after stem analysis is completed.
#'
#' @return Result of stem analysis stored in a tibble reporting all possible
#'   model profiles with p-values, associated gene members and more...
run_stem <- function(path, jar_path = "external_software/stem/stem.jar",
                     clear_output = F) {
  input_path <- file.path(path, "input")
  setting_path <- file.path(str_replace(path, "output", "data"), "setting")
  out_path <- file.path(path, "output")

  # stem translates all gene symbols to upper case so we need to create a
  # mapping between the original gene id and the upper case version
  gene_mapping <- list.files(input_path, full.names = T, pattern = ".txt") %>%
    map_dfr(function(p) {
      read_delim(p, delim = "\t") %>%
        distinct(symbol = gene)
    }) %>%
    mutate(gene = str_to_upper(symbol)) %>%
    distinct()

  stem_call <- str_c("java -mx1024M -jar", jar_path, "-b", setting_path,
    out_path,
    sep = " "
  )
  # run stem
  system(stem_call)

  # combine stem results
  # get gene tables
  gtables <- list.files(out_path, pattern = "genetable", full.names = T) %>%
    map_dfr(function(g) {
      key <- basename(g) %>%
        str_remove("_genetable.txt")
      read_delim(g, delim = "\t") %>%
        mutate(key = key) %>%
        janitor::clean_names() %>%
        select(-spot) %>%
        gather(time, value, -gene, -profile, -key) %>%
        mutate(time = parse_number(time))
    })

  # get profile tables
  ptables <- list.files(out_path,
    pattern = "profiletable",
    full.names = T
  ) %>%
    map_dfr(function(g) {
      key <- basename(g) %>%
        str_remove("_profiletable.txt")
      read_delim(g, delim = "\t") %>%
        mutate(key = key) %>%
        janitor::clean_names() %>%
        select(
          profile = profile_id, y_coords = profile_model,
          size = number_genes_assigned, p = p_value, key
        )
    })

  # combine tables
  tables <- inner_join(gtables, ptables) %>%
    inner_join(gene_mapping, by = "gene") %>%
    transmute(key, profile, size, p, y_coords, gene = symbol, time, value)

  if (clear_output) {
    # delete output files from stem
    file.remove(list.files(out_path, full.names = T))
  }
  return(tables)
}

#' This functions loads gene sets from GO-BP, progeny and dorothea for a
#' specified organism
#'
#' @param organism character; either "mouse" or "human".
#'
#' @return Tibble of gene sets.
load_genesets <- function(organism = "mouse") {
  if (organism == "mouse") {
    go_terms <- msigdf::msigdf.mouse %>%
      filter(category_subcode == "bp") %>%
      select(geneset, gene = mouse.symbol) %>%
      mutate(group = "go")

    progeny <- progeny::getModel(organism = "Mouse", top = 100) %>%
      tdy("gene", "geneset", "weight") %>%
      filter(weight != 0) %>%
      select(gene, geneset) %>%
      mutate(group = "progeny")

    dorothea <- dorothea::dorothea_mm %>%
      select(geneset = tf, gene = target, confidence) %>%
      mutate(group = "dorothea")
  } else if (organism == "human") {
    go_terms <- msigdf::msigdf.human %>%
      filter(category_subcode == "bp") %>%
      select(geneset, gene = symbol) %>%
      mutate(group = "go")

    progeny <- progeny::getModel(organism = "Human", top = 100) %>%
      tdy("gene", "geneset", "weight") %>%
      filter(weight != 0) %>%
      select(gene, geneset) %>%
      mutate(group = "progeny")

    dorothea <- dorothea::dorothea_hs %>%
      select(geneset = tf, gene = target, confidence) %>%
      mutate(group = "dorothea")
  }

  genesets <- bind_rows(go_terms, progeny, dorothea)

  return(genesets)
}
