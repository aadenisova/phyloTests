#!/usr/bin/env Rscript
# Aggregate per-phenotype result tables into master summaries.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

PROJECT_ROOT <- normalizePath(".", wins = FALSE)
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
MASTER_DIR <- file.path(RESULTS_DIR, "_master")

PHENOTYPE_MAP <- c(
  TOTAL = "TOTALINNOVATIONS2025_ResEff",
  FOOD = "FOODINNO2025_ResEff",
  TECH = "TECHINNO2025_ResEff"
)
INNOVATION_TRAITS <- unname(PHENOTYPE_MAP)

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

read_phenotype_table <- function(phenotype, filename) {
  path <- file.path(RESULTS_DIR, PHENOTYPE_MAP[[phenotype]], filename)
  if (!file.exists(path)) {
    stop("Missing file: ", path)
  }
  df <- read.delim(path, stringsAsFactors = FALSE)
  df$phenotype <- phenotype
  df
}

pivot_beta_matrix <- function(df, value_col = "beta") {
  df |>
    select(trait, phenotype, value = all_of(value_col)) |>
    pivot_wider(names_from = phenotype, values_from = value) |>
    arrange(trait)
}

pivot_lambda_matrix <- function(df) {
  innov <- df |>
    filter(trait %in% INNOVATION_TRAITS) |>
    select(trait, phenotype, lambda)

  innov |>
    mutate(innovation = recode(
      trait,
      TOTALINNOVATIONS2025_ResEff = "TOTAL",
      FOODINNO2025_ResEff = "FOOD",
      TECHINNO2025_ResEff = "TECH"
    )) |>
    select(innovation, phenotype, lambda) |>
    pivot_wider(names_from = phenotype, values_from = lambda) |>
    arrange(innovation)
}

# ---- master tables ----
phylo_master <- bind_rows(lapply(names(PHENOTYPE_MAP), function(ph) {
  read_phenotype_table(ph, "phylogenetic_signal.tsv") |>
    transmute(
      trait = Trait,
      K, K_p, lambda, lambda_p,
      phenotype
    )
}))

pgls_uni_master <- bind_rows(lapply(names(PHENOTYPE_MAP), function(ph) {
  read_phenotype_table(ph, "pgls_univariate.tsv") |>
    transmute(
      trait = Trait,
      beta, se = SE, p, r2 = R2,
      phenotype
    )
}))

pgls_mv_master <- bind_rows(lapply(names(PHENOTYPE_MAP), function(ph) {
  read_phenotype_table(ph, "pgls_multivariate.tsv") |>
    transmute(
      trait = term,
      beta, se = SE, p, r2 = model_R2,
      phenotype
    )
}))

top20_master <- bind_rows(lapply(names(PHENOTYPE_MAP), function(ph) {
  read_phenotype_table(ph, "top20_branch_changes.tsv")
}))

# ---- summary matrices ----
pgls_uni_beta_matrix <- pivot_beta_matrix(pgls_uni_master)
pgls_mv_beta_matrix <- pivot_beta_matrix(pgls_mv_master)
innovation_lambda_matrix <- pivot_lambda_matrix(phylo_master)

# ---- save ----
dir.create(MASTER_DIR, recursive = TRUE, showWarnings = FALSE)

write_tsv(phylo_master, file.path(MASTER_DIR, "phylogenetic_signal_master.tsv"))
write_tsv(pgls_uni_master, file.path(MASTER_DIR, "pgls_univariate_master.tsv"))
write_tsv(pgls_mv_master, file.path(MASTER_DIR, "pgls_multivariate_master.tsv"))
write_tsv(top20_master, file.path(MASTER_DIR, "top20_branch_changes_master.tsv"))
write_tsv(pgls_uni_beta_matrix, file.path(MASTER_DIR, "pgls_univariate_beta_matrix.tsv"))
write_tsv(pgls_mv_beta_matrix, file.path(MASTER_DIR, "pgls_multivariate_beta_matrix.tsv"))
write_tsv(innovation_lambda_matrix, file.path(MASTER_DIR, "innovation_lambda_matrix.tsv"))

cat("Master tables written to:", MASTER_DIR, "\n")
cat("  phylogenetic_signal:", nrow(phylo_master), "rows\n")
cat("  pgls_univariate:", nrow(pgls_uni_master), "rows\n")
cat("  pgls_multivariate:", nrow(pgls_mv_master), "rows\n")
cat("  top20_branch_changes:", nrow(top20_master), "rows\n")
