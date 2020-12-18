library(biomaRt)
library(tidyverse)
library(here)
library(GO.db)

# mapping of gene symbols -------------------------------------------------
# download gene annotations across the species homo sapiens and mus musculus for
# various gene identifiers (e.g. symbol, ensembl, ...)

# listEnsemblArchives()
# used version/host = "http://jan2020.archive.ensembl.org"

mouse_ensembl <- useMart("ensembl",
  dataset = "mmusculus_gene_ensembl",
  host = "http://jan2020.archive.ensembl.org"
)
human_ensembl <- useMart("ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "http://jan2020.archive.ensembl.org"
)

common_attributes <- c(
  "ensembl_gene_id", "ensembl_gene_id_version",
  "entrezgene_id"
)

biomart_output <- getLDS(
  attributes = c("mgi_symbol", common_attributes), mart = mouse_ensembl,
  attributesL = c("hgnc_symbol", common_attributes),
  martL = human_ensembl
)

tidy_biomart_output <- biomart_output %>%
  as_tibble() %>%
  rename(
    symbol_mgi = MGI.symbol,
    ensembl_mgi = Gene.stable.ID,
    ensembl_v_mgi = Gene.stable.ID.version,
    entrez_mgi = NCBI.gene.ID,
    symbol_hgnc = HGNC.symbol,
    ensembl_hgnc = Gene.stable.ID.1,
    ensembl_v_hgnc = Gene.stable.ID.version.1,
    entrez_hgnc = NCBI.gene.ID.1
  ) %>%
  mutate_if(is.integer, as.character) %>%
  na_if("")

saveRDS(
  tidy_biomart_output,
  here("data/annotation/gene_id_annotation.rds")
)

# microarray platforms ----------------------------------------------------
# Mapping of annotation identifier to annotation package
#
# there should be a mapping between the two available somewhere?
# List of BioC annots:
# https://www.bioconductor.org/packages/release/data/annotation/
mapping <- list(
  "pd.mouse430.2" = "mouse4302.db",
  "pd.mogene.2.0.st" = "mogene20sttranscriptcluster.db",
  "pd.hg.u133.plus.2" = "hgu133plus2.db",
  "pd.hg.u219" = "hgu219.db",
  "pd.hugene.1.1.st.v1" = "hugene11sttranscriptcluster.db",
  "pd.hg.u133a" = "hgu133a.db"
)

saveRDS(mapping, here("data/annotation/platforms.rds"))


# go term mapping ---------------------------------------------------------
go_mapping <- Term(GOTERM) %>%
  enframe("id", "term") %>%
  mutate(
    clean_term = str_to_lower(term),
    clean_term = str_replace_all(clean_term, "-", " "),
    clean_term = str_replace_all(clean_term, "/", " "),
    clean_term = str_replace_all(clean_term, "\\(", ""),
    clean_term = str_replace_all(clean_term, "\\)", ""),
    clean_term = str_replace_all(clean_term, ",", ""),
    clean_term = str_replace_all(clean_term, "'", "")
  )

saveRDS(go_mapping, here("data/annotation/go_mapping.rds"))

# go stop words -----------------------------------------------------------
# word to ignore for semantic analysis of GO terms
stop_go_words <- tibble(
  word = c(
    "cell", "cellular", "positive", "negative", "periphery", "mediated",
    "regulation", "process", "establishment", "involved", "particle",
    "complex", "substance", "compound", "system", "response", "protein",
    "beta", "induced", "gamma", "incorrect", "class", "processing",
    "family", "chain", "alpha", "activity", "stimulus", "projection",
    "pathway", "component", "activation", "phase", "transition",
    "production", "organization", "biological", "based", "structure",
    "extension", "levels", "modified"
  )
) %>%
  distinct(word)

saveRDS(stop_go_words, here("data/annotation/stop_go_words.rds"))
