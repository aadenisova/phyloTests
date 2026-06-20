#!/usr/bin/env Rscript
# Phylogenetic analysis of bird innovativeness (phenotype + tree only).

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(phylolm)
  library(ggplot2)
  library(ggtree)
  library(dplyr)
  library(tidyr)
})

# ---- parameters ----
PROJECT_ROOT <- normalizePath(".", wins = FALSE)
DATA_DIR <- file.path(PROJECT_ROOT, "data")
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
TREE_FILE <- "roadies_birds_allbirdtraits.nwk"
TRAITS_FILE <- "ALLBIRDTRAITS_intersect.csv"
SPECIES_COL <- "Species"

PHENOTYPES <- c(
  "TOTALINNOVATIONS2025_ResEff",
  "FOODINNO2025_ResEff",
  "TECHINNO2025_ResEff"
)
COVARIATES <- c(
  "Brain_size", "Relative_brain_size", "Mass", "DietBreadth",
  "HabitatBreadth", "GenerationLength", "Range.Size",
  "UrbanFULL", "Migration", "Trophic_level"
)
SIGNAL_TRAITS <- c(
  PHENOTYPES,
  "Brain_size", "Relative_brain_size", "Mass",
  "DietBreadth", "HabitatBreadth", "GenerationLength", "Range.Size"
)
MULTIVARIATE_COVARIATES <- c(
  "Relative_brain_size", "Mass", "DietBreadth",
  "HabitatBreadth", "GenerationLength"
)
PCA_TRAITS <- c(
  "Relative_brain_size", "Mass", "DietBreadth",
  "HabitatBreadth", "GenerationLength", "Range.Size"
)
CORRELATION_TRAITS <- c(
  PHENOTYPES,
  "Brain_size", "Relative_brain_size", "Mass",
  "DietBreadth", "HabitatBreadth", "GenerationLength", "Range.Size"
)
FOREGROUND_SD <- 1
# PGLS predictors are z-scaled inside pgls_fit() for numerical stability
# (Range.Size otherwise causes a singular design matrix in phylolm).

# ---- helpers ----
dir_ok <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

load_inputs <- function() {
  tree <- read.tree(file.path(DATA_DIR, TREE_FILE))
  dat <- read.csv(
    file.path(DATA_DIR, TRAITS_FILE),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if ("Unnamed: 0" %in% names(dat)) {
    dat <- dat[, setdiff(names(dat), "Unnamed: 0"), drop = FALSE]
  }
  list(tree = tree, dat = dat)
}

prepare_data <- function(tree, dat) {
  tree_species <- sort(tree$tip.label)
  table_species <- sort(unique(dat[[SPECIES_COL]]))

  in_both <- intersect(tree_species, table_species)
  missing_in_table <- setdiff(tree_species, table_species)
  missing_in_tree <- setdiff(table_species, tree_species)

  dat <- dat[dat[[SPECIES_COL]] %in% in_both, , drop = FALSE]
  dat <- dat[match(tree$tip.label, dat[[SPECIES_COL]]), , drop = FALSE]
  rownames(dat) <- dat[[SPECIES_COL]]

  report <- data.frame(
    metric = c(
      "species_in_tree", "species_in_table",
      "species_matched", "missing_in_table", "missing_in_tree"
    ),
    value = c(
      length(tree_species), length(table_species),
      length(in_both), length(missing_in_table), length(missing_in_tree)
    ),
    stringsAsFactors = FALSE
  )

  list(
    tree = tree,
    dat = dat,
    report = report,
    missing_in_table = missing_in_table,
    missing_in_tree = missing_in_tree
  )
}

trait_vector <- function(dat, tree, trait) {
  x <- dat[[trait]]
  names(x) <- dat[[SPECIES_COL]]
  x <- x[tree$tip.label]
  x
}

prune_complete <- function(tree, dat, traits) {
  ok <- complete.cases(dat[, traits, drop = FALSE])
  tree_pr <- drop.tip(tree, tree$tip.label[!ok])
  dat_pr <- dat[tree_pr$tip.label, , drop = FALSE]
  list(tree = tree_pr, dat = dat_pr)
}

compute_phylogenetic_signal <- function(tree, dat, traits) {
  rows <- lapply(traits, function(tr) {
    cc <- prune_complete(tree, dat, tr)
    x <- trait_vector(cc$dat, cc$tree, tr)
    k <- phytools::phylosig(cc$tree, x, method = "K", test = TRUE)
    l <- phytools::phylosig(cc$tree, x, method = "lambda", test = TRUE)
    data.frame(
      Trait = tr,
      K = unname(k$K),
      K_p = k$P,
      lambda = unname(l$lambda),
      lambda_p = l$P,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

compute_ancestral_states <- function(tree, x) {
  anc <- fastAnc(tree, x)
  data.frame(
    node = as.integer(names(anc)),
    state = as.numeric(anc),
    stringsAsFactors = FALSE
  )
}

compute_branch_changes <- function(tree, x, anc) {
  ntip <- Ntip(tree)
  edge <- tree$edge
  tip_vals <- setNames(as.numeric(x), names(x))

  node_value <- function(node) {
    if (node <= ntip) {
      return(tip_vals[tree$tip.label[node]])
    }
    anc[[as.character(node)]]
  }

  out <- data.frame(
    branch_id = paste0("edge_", seq_len(nrow(edge))),
    parent = edge[, 1],
    child = edge[, 2],
    parent_value = vapply(edge[, 1], node_value, numeric(1)),
    child_value = vapply(edge[, 2], node_value, numeric(1)),
    branch_length = tree$edge.length,
    stringsAsFactors = FALSE
  )
  out$delta <- out$child_value - out$parent_value
  out
}

identify_foreground <- function(branch_changes, sd_mult = FOREGROUND_SD) {
  threshold <- mean(branch_changes$delta) + sd_mult * sd(branch_changes$delta)
  fg <- branch_changes[branch_changes$delta > threshold, , drop = FALSE]
  list(foreground = fg, threshold = threshold)
}

pgls_fit <- function(tree, dat, formula, scale_terms = NULL) {
  d <- dat
  if (!is.null(scale_terms)) {
    for (tr in scale_terms) {
      if (tr %in% names(d) && is.numeric(d[[tr]])) {
        d[[tr]] <- as.numeric(scale(d[[tr]]))
      }
    }
  }
  phylolm(
    formula,
    data = d,
    phy = tree,
    model = "BM",
    boot = 0
  )
}

extract_slope <- function(fit, term) {
  coefs <- summary(fit)$coefficients
  if (!term %in% rownames(coefs)) {
    return(data.frame(beta = NA_real_, SE = NA_real_, p = NA_real_, R2 = NA_real_))
  }
  data.frame(
    beta = coefs[term, "Estimate"],
    SE = coefs[term, "StdErr"],
    p = coefs[term, "p.value"],
    R2 = summary(fit)$adj.r.squared,
    stringsAsFactors = FALSE
  )
}

run_univariate_pgls <- function(tree, dat, phenotype, covariates) {
  rows <- lapply(covariates, function(tr) {
    d <- dat
    if (tr %in% c("Migration", "Trophic_level")) {
      d[[tr]] <- as.numeric(factor(d[[tr]]))
    }
    fit <- tryCatch(
      pgls_fit(tree, d, as.formula(paste(phenotype, "~", tr)), scale_terms = tr),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.frame(Trait = tr, beta = NA_real_, SE = NA_real_, p = NA_real_, R2 = NA_real_))
    }
    out <- extract_slope(fit, tr)
    out$Trait <- tr
    out
  })
  bind_rows(rows) |>
    transmute(Trait, beta, SE, p, R2) |>
    arrange(p)
}

run_multivariate_pgls <- function(tree, dat, phenotype, predictors) {
  rhs <- paste(predictors, collapse = " + ")
  fit <- tryCatch(
    pgls_fit(
      tree, dat,
      as.formula(paste(phenotype, "~", rhs)),
      scale_terms = predictors
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    out <- data.frame(
      term = predictors,
      beta = NA_real_,
      SE = NA_real_,
      p = NA_real_,
      stringsAsFactors = FALSE
    )
    attr(out, "R2") <- NA_real_
    return(out)
  }
  coefs <- summary(fit)$coefficients
  keep <- rownames(coefs) != "(Intercept)"
  out <- data.frame(
    term = rownames(coefs)[keep],
    beta = coefs[keep, "Estimate"],
    SE = coefs[keep, "StdErr"],
    p = coefs[keep, "p.value"],
    stringsAsFactors = FALSE
  )
  attr(out, "R2") <- summary(fit)$adj.r.squared
  out
}

run_pca <- function(dat, traits, phenotype) {
  cc <- complete.cases(dat[, traits, drop = FALSE])
  x <- scale(dat[cc, traits, drop = FALSE])
  pca <- prcomp(x, center = FALSE, scale. = FALSE)
  coords <- data.frame(
    species = dat[[SPECIES_COL]][cc],
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    innovation = dat[[phenotype]][cc],
    stringsAsFactors = FALSE
  )
  list(coords = coords, pca = pca)
}

save_genome_inputs <- function(dat, phenotype, out_dir) {
  write_tsv(
    data.frame(
      species = dat[[SPECIES_COL]],
      innovation = dat[[phenotype]],
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "continuous_phenotype.tsv")
  )
}

# ---- figures ----
plot_contmap <- function(tree, x, out_file) {
  pdf(out_file, width = 12, height = 16)
  contMap(
    tree, x,
    lut = c("blue", "white", "red"),
    outline = FALSE,
    lwd = 3,
    fsize = 0.35
  )
  dev.off()
}

plot_ancestral_tree <- function(tree, anc_df, out_file) {
  anc_df$label <- round(anc_df$state, 2)
  p <- ggtree(tree, branch.length = "none", layout = "rectangular") %<+% anc_df +
    geom_nodelab(aes(label = label), size = 1.2, color = "darkblue") +
    theme_tree2() +
    ggtitle("Ancestral state reconstruction")
  ggsave(out_file, p, width = 12, height = 16)
}

plot_foreground_tree <- function(tree, branch_changes, foreground, out_file) {
  ne <- Nedge(tree)
  edge_col <- rep("gray60", ne)
  edge_w <- rep(0.6, ne)
  fg_ids <- foreground$branch_id
  idx <- match(fg_ids, branch_changes$branch_id)
  idx <- idx[!is.na(idx)]
  edge_col[idx] <- "red"
  edge_w[idx] <- 2.5

  pdf(out_file, width = 12, height = 16)
  plot(
    tree,
    type = "phylogram",
    edge.color = edge_col,
    edge.width = edge_w,
    cex = 0.35,
    no.margin = TRUE
  )
  title(main = "Foreground branches (independent innovation candidates)")
  dev.off()
}

plot_delta_histogram <- function(branch_changes, threshold, out_file) {
  p <- ggplot(branch_changes, aes(delta)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "white") +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = "Distribution of branch-wise trait changes",
      x = expression(Delta ~ "(child - parent)"),
      y = "Count"
    ) +
    theme_bw()
  ggsave(out_file, p, width = 7, height = 5)
}

plot_parent_child_scatter <- function(branch_changes, foreground, out_file) {
  branch_changes$foreground <- branch_changes$branch_id %in% foreground$branch_id
  p <- ggplot(branch_changes, aes(parent_value, child_value, color = foreground)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red")) +
    labs(
      title = "Parent vs child values along branches",
      x = "Parent value",
      y = "Child value",
      color = "Foreground"
    ) +
    theme_bw()
  ggsave(out_file, p, width = 6, height = 5)
}

plot_forest <- function(pgls_uni, out_file) {
  pgls_uni <- pgls_uni |>
    mutate(
      lower = beta - 1.96 * SE,
      upper = beta + 1.96 * SE
    )
  p <- ggplot(pgls_uni, aes(x = beta, y = reorder(Trait, beta))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.2) +
    geom_point(size = 2) +
    labs(
      title = "Univariate PGLS coefficients (95% CI)",
      x = expression(beta),
      y = NULL
    ) +
    theme_bw()
  ggsave(out_file, p, width = 7, height = 5)
}

plot_correlation_matrix <- function(dat, traits, out_file) {
  cor_mat <- cor(dat[, traits, drop = FALSE], use = "pairwise.complete.obs")
  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c("Var1", "Var2", "r")
  p <- ggplot(cor_df, aes(Var1, Var2, fill = r)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    coord_fixed() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7)
    ) +
    labs(title = "Trait correlation matrix", x = NULL, y = NULL, fill = "r")
  ggsave(out_file, p, width = 9, height = 8)
}

plot_pca <- function(pca_coords, out_file) {
  p <- ggplot(pca_coords, aes(PC1, PC2, color = innovation)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = median(pca_coords$innovation, na.rm = TRUE)
    ) +
    labs(
      title = "PCA of ecological traits",
      x = "PC1",
      y = "PC2",
      color = "Innovation"
    ) +
    theme_bw()
  ggsave(out_file, p, width = 7, height = 6)
}

# ---- per-phenotype pipeline ----
run_phenotype_analysis <- function(tree, dat, phenotype, prep_info) {
  out_dir <- file.path(RESULTS_DIR, phenotype)
  fig_dir <- file.path(out_dir, "figures")
  dir_ok(out_dir)
  dir_ok(fig_dir)

  message("=== ", phenotype, " ===")

  # Step 1
  write_tsv(dat, file.path(out_dir, "cleaned_data.tsv"))
  write_tsv(prep_info$report, file.path(out_dir, "data_match_report.tsv"))
  write_tsv(
    data.frame(species = prep_info$missing_in_table, stringsAsFactors = FALSE),
    file.path(out_dir, "missing_in_table.tsv")
  )
  write_tsv(
    data.frame(species = prep_info$missing_in_tree, stringsAsFactors = FALSE),
    file.path(out_dir, "missing_in_tree.tsv")
  )
  cat(
    phenotype, ": tree =", prep_info$report$value[1],
    ", table =", prep_info$report$value[2],
    ", matched =", prep_info$report$value[3], "\n"
  )

  x <- trait_vector(dat, tree, phenotype)

  # Step 2
  signal <- compute_phylogenetic_signal(tree, dat, SIGNAL_TRAITS)
  write_tsv(signal, file.path(out_dir, "phylogenetic_signal.tsv"))

  # Step 3
  anc_df <- compute_ancestral_states(tree, x)
  write_tsv(anc_df, file.path(out_dir, "ancestral_states.tsv"))
  anc_named <- setNames(anc_df$state, anc_df$node)

  # Step 4
  branch_changes <- compute_branch_changes(tree, x, anc_named)
  write_tsv(branch_changes, file.path(out_dir, "branch_changes.tsv"))

  # Step 5
  fg_info <- identify_foreground(branch_changes, FOREGROUND_SD)
  write_tsv(fg_info$foreground, file.path(out_dir, "foreground_branches.tsv"))
  top20 <- branch_changes |>
    arrange(desc(delta)) |>
    slice_head(n = min(20, nrow(branch_changes)))
  write_tsv(top20, file.path(out_dir, "top20_branch_changes.tsv"))

  # Step 6
  pgls_uni <- run_univariate_pgls(tree, dat, phenotype, COVARIATES)
  write_tsv(pgls_uni, file.path(out_dir, "pgls_univariate.tsv"))

  # Step 7
  cc <- complete.cases(dat[, c(phenotype, MULTIVARIATE_COVARIATES), drop = FALSE])
  tree_mv <- drop.tip(tree, tree$tip.label[!cc])
  dat_mv <- dat[tree_mv$tip.label, , drop = FALSE]
  pgls_mv <- run_multivariate_pgls(tree_mv, dat_mv, phenotype, MULTIVARIATE_COVARIATES)
  pgls_mv_out <- transform(
    pgls_mv,
    model_R2 = attr(pgls_mv, "R2")
  )
  write_tsv(pgls_mv_out, file.path(out_dir, "pgls_multivariate.tsv"))

  # Step 8
  pca_res <- run_pca(dat, PCA_TRAITS, phenotype)
  write_tsv(pca_res$coords, file.path(out_dir, "pca_coordinates.tsv"))

  # Step 9
  save_genome_inputs(dat, phenotype, out_dir)

  # Figures
  plot_contmap(tree, x, file.path(fig_dir, "contMap.pdf"))
  plot_ancestral_tree(tree, anc_df, file.path(fig_dir, "ancestral_states.pdf"))
  plot_foreground_tree(
    tree, branch_changes, fg_info$foreground,
    file.path(fig_dir, "foreground_branches.pdf")
  )
  plot_delta_histogram(
    branch_changes, fg_info$threshold,
    file.path(fig_dir, "delta_histogram.pdf")
  )
  plot_parent_child_scatter(
    branch_changes, fg_info$foreground,
    file.path(fig_dir, "parent_child_scatter.pdf")
  )
  plot_forest(pgls_uni, file.path(fig_dir, "forest_plot.pdf"))
  plot_correlation_matrix(
    dat, CORRELATION_TRAITS,
    file.path(fig_dir, "correlation_matrix.pdf")
  )
  plot_pca(pca_res$coords, file.path(fig_dir, "pca.pdf"))

  invisible(out_dir)
}

# ---- main ----
main <- function() {
  inputs <- load_inputs()
  prep <- prepare_data(inputs$tree, inputs$dat)

  lapply(PHENOTYPES, function(ph) {
    run_phenotype_analysis(prep$tree, prep$dat, ph, prep)
  })

  message("Done. Results in: ", RESULTS_DIR)
}

if (sys.nframe() == 0) {
  main()
}
