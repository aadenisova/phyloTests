#!/usr/bin/env Rscript
# Exploratory trait–innovation correlations on the PGLS-ready subset.
# Reads precomputed cleaned_data.tsv; does not modify the main pipeline.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---- parameters ----
PROJECT_ROOT <- normalizePath(".", wins = FALSE)
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
OUT_DIR <- file.path(RESULTS_DIR, "explore_traits")
FIG_DIR <- file.path(OUT_DIR, "figures")

# Same table used by run_phenotype_analysis() (83 matched tree tips).
CLEANED_DATA <- file.path(
  RESULTS_DIR, "TOTALINNOVATIONS2025_ResEff", "cleaned_data.tsv"
)

PHENOTYPES <- c(
  "TOTALINNOVATIONS2025",
  "TOTALINNOVATIONS2025_ResEff",
  "FOODINNO2025",
  "FOODINNO2025_ResEff",
  "TECHINNO2025",
  "TECHINNO2025_ResEff"
)

PREDICTORS <- c(
  "Brain_size", "Relative_brain_size", "Mass", "GenerationLength",
  "Range.Size", "HabitatBreadth", "DietBreadth",
  "UrbanFULL", "Migration", "Trophic_level"
)

CONTINUOUS_PREDICTORS <- c(
  "Brain_size", "Relative_brain_size", "Mass", "GenerationLength",
  "Range.Size", "HabitatBreadth"
)

ORDINAL_PREDICTORS <- "DietBreadth"

CATEGORICAL_PREDICTORS <- c("UrbanFULL", "Migration", "Trophic_level")

CORRELATION_MATRIX_TRAITS <- c(
  PHENOTYPES,
  "Brain_size", "Relative_brain_size", "Mass",
  "GenerationLength", "Range.Size", "HabitatBreadth"
)

SIGN_FLIP_THRESHOLD <- 0.05
NONLINEARITY_THRESHOLD <- 0.15
OUTLIER_STD_RESID <- 2.5

# ---- helpers ----
dir_ok <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", row.names = FALSE, quote = FALSE)
}

load_cleaned_data <- function(path = CLEANED_DATA) {
  if (!file.exists(path)) {
    stop("cleaned_data.tsv not found: ", path)
  }
  dat <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  for (col in c(PHENOTYPES, PREDICTORS, CORRELATION_MATRIX_TRAITS)) {
    if (col %in% names(dat) && !is.numeric(dat[[col]])) {
      dat[[col]] <- suppressWarnings(as.numeric(dat[[col]]))
    }
  }
  message("Loaded ", nrow(dat), " species from ", path)
  dat
}

predictor_numeric <- function(x, trait) {
  if (trait %in% c("Migration", "Trophic_level")) {
    return(as.numeric(factor(x)))
  }
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  suppressWarnings(as.numeric(x))
}

pairwise_complete <- function(y, x) {
  ok <- complete.cases(y, x)
  list(y = y[ok], x = x[ok], n = sum(ok))
}

compute_correlation_row <- function(phenotype, predictor, dat) {
  y <- dat[[phenotype]]
  x <- predictor_numeric(dat[[predictor]], predictor)
  cc <- pairwise_complete(y, x)
  if (cc$n < 3) {
    return(data.frame(
      phenotype = phenotype,
      predictor = predictor,
      pearson_r = NA_real_,
      spearman_rho = NA_real_,
      pearson_p = NA_real_,
      spearman_p = NA_real_,
      n = cc$n,
      stringsAsFactors = FALSE
    ))
  }
  pearson <- suppressWarnings(cor.test(cc$y, cc$x, method = "pearson"))
  spearman <- suppressWarnings(cor.test(cc$y, cc$x, method = "spearman"))
  data.frame(
    phenotype = phenotype,
    predictor = predictor,
    pearson_r = unname(pearson$estimate),
    spearman_rho = unname(spearman$estimate),
    pearson_p = pearson$p.value,
    spearman_p = spearman$p.value,
    n = cc$n,
    stringsAsFactors = FALSE
  )
}

compute_all_correlations <- function(dat) {
  rows <- lapply(PHENOTYPES, function(ph) {
    bind_rows(lapply(PREDICTORS, function(pr) {
      compute_correlation_row(ph, pr, dat)
    }))
  })
  bind_rows(rows)
}

format_stats_label <- function(pearson_r, spearman_rho, pearson_p) {
  sprintf(
    "Pearson r = %.3f\nSpearman rho = %.3f\np = %.3g",
    pearson_r, spearman_rho, pearson_p
  )
}

plot_continuous_scatter <- function(dat, phenotype, predictor, stats_row, out_file) {
  df <- data.frame(
    innovation = dat[[phenotype]],
    trait = dat[[predictor]],
    species = dat[["Species"]],
    stringsAsFactors = FALSE
  )
  df <- df[complete.cases(df[, c("innovation", "trait")]), , drop = FALSE]

  p <- ggplot(df, aes(trait, innovation)) +
    geom_point(size = 2, alpha = 0.75, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "firebrick", fill = "pink") +
    labs(
      title = paste0(phenotype, " vs ", predictor),
      x = predictor,
      y = phenotype
    ) +
    annotate(
      "label",
      x = Inf, y = Inf,
      label = format_stats_label(
        stats_row$pearson_r, stats_row$spearman_rho, stats_row$pearson_p
      ),
      hjust = 1.05, vjust = 1.1,
      size = 3.2,
      label.size = 0.2,
      fill = "white",
      alpha = 0.85
    ) +
    theme_bw()
  ggsave(out_file, p, width = 6, height = 5)
}

plot_dietbreadth <- function(dat, phenotype, out_dir) {
  df <- data.frame(
    innovation = dat[[phenotype]],
    diet = dat[["DietBreadth"]],
    stringsAsFactors = FALSE
  )
  df <- df[complete.cases(df), , drop = FALSE]
  df$diet_f <- factor(df$diet)

  p_box <- ggplot(df, aes(diet_f, innovation)) +
    geom_boxplot(fill = "lightblue", color = "gray30", outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = 0.7, size = 1.8, color = "steelblue") +
    labs(
      title = paste0(phenotype, " by DietBreadth"),
      x = "DietBreadth",
      y = phenotype
    ) +
    theme_bw()
  ggsave(
    file.path(out_dir, paste0("DietBreadth_boxplot_", phenotype, ".pdf")),
    p_box, width = 6, height = 5
  )

  p_scatter <- ggplot(df, aes(diet, innovation)) +
    geom_point(
      position = position_jitter(width = 0.15, height = 0),
      size = 2, alpha = 0.75, color = "steelblue"
    ) +
    geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "firebrick", fill = "pink") +
    labs(
      title = paste0(phenotype, " vs DietBreadth (jittered)"),
      x = "DietBreadth",
      y = phenotype
    ) +
    theme_bw()
  ggsave(
    file.path(out_dir, paste0("DietBreadth_scatter_", phenotype, ".pdf")),
    p_scatter, width = 6, height = 5
  )
}

plot_categorical <- function(dat, phenotype, predictor, out_dir) {
  raw <- dat[[predictor]]
  df <- data.frame(
    innovation = dat[[phenotype]],
    group = raw,
    stringsAsFactors = FALSE
  )
  df <- df[complete.cases(df$innovation) & !is.na(df$group) & df$group != "", , drop = FALSE]
  df$group <- factor(df$group)

  p_box <- ggplot(df, aes(group, innovation)) +
    geom_boxplot(fill = "lightblue", color = "gray30", outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.8, color = "steelblue") +
    labs(
      title = paste0(phenotype, " by ", predictor),
      x = predictor,
      y = phenotype
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(
    file.path(out_dir, paste0(predictor, "_boxplot_", phenotype, ".pdf")),
    p_box, width = 6.5, height = 5
  )

  p_violin <- ggplot(df, aes(group, innovation)) +
    geom_violin(fill = "lightblue", color = "gray30", alpha = 0.65) +
    geom_jitter(width = 0.12, alpha = 0.75, size = 1.8, color = "steelblue") +
    labs(
      title = paste0(phenotype, " by ", predictor, " (violin)"),
      x = predictor,
      y = phenotype
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(
    file.path(out_dir, paste0(predictor, "_violin_", phenotype, ".pdf")),
    p_violin, width = 6.5, height = 5
  )
}

plot_correlation_matrix <- function(dat, traits, out_file) {
  mat <- cor(dat[, traits, drop = FALSE], use = "pairwise.complete.obs")
  cor_df <- as.data.frame(as.table(mat))
  names(cor_df) <- c("Var1", "Var2", "r")
  cor_df$Var1 <- factor(cor_df$Var1, levels = traits)
  cor_df$Var2 <- factor(cor_df$Var2, levels = rev(traits))

  p <- ggplot(cor_df, aes(Var1, Var2, fill = r)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", r)), size = 2.2, color = "black") +
    scale_fill_gradient2(
      low = "#2166ac", mid = "white", high = "#b2182b",
      midpoint = 0, limits = c(-1, 1), name = "r"
    ) +
    coord_fixed() +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 7)
    ) +
    labs(
      title = "Innovation and trait correlation matrix (Pearson, pairwise complete)",
      x = NULL, y = NULL
    )
  ggsave(out_file, p, width = 11, height = 9)
}

detect_outliers <- function(dat, phenotype, predictor) {
  y <- dat[[phenotype]]
  x <- predictor_numeric(dat[[predictor]], predictor)
  ok <- complete.cases(y, x)
  if (sum(ok) < 4) {
    return(character(0))
  }
  fit <- lm(y[ok] ~ x[ok])
  std_resid <- rstandard(fit)
  idx <- which(abs(std_resid) > OUTLIER_STD_RESID)
  if (length(idx) == 0) {
    return(character(0))
  }
  dat[["Species"]][ok][idx]
}

detect_nonlinearity <- function(dat, phenotype, predictor) {
  if (!predictor %in% CONTINUOUS_PREDICTORS) {
    return(FALSE)
  }
  y <- dat[[phenotype]]
  x <- dat[[predictor]]
  cc <- pairwise_complete(y, x)
  if (cc$n < 6) {
    return(FALSE)
  }
  lin <- lm(cc$y ~ cc$x)
  quad <- lm(cc$y ~ cc$x + I(cc$x^2))
  lin_r2 <- summary(lin)$r.squared
  quad_r2 <- summary(quad)$r.squared
  (quad_r2 - lin_r2) > 0.08 && anova(lin, quad)$`Pr(>F)`[2] < 0.1
}

format_top_list <- function(rows, col, direction = c("positive", "negative"), n = 3,
                            p_col = "pearson_p") {
  direction <- match.arg(direction)
  rows <- rows[!is.na(rows[[col]]), , drop = FALSE]
  if (direction == "positive") {
    rows <- rows[order(-rows[[col]], rows[[p_col]]), , drop = FALSE]
    rows <- rows[rows[[col]] > 0, , drop = FALSE]
  } else {
    rows <- rows[order(rows[[col]], rows[[p_col]]), , drop = FALSE]
    rows <- rows[rows[[col]] < 0, , drop = FALSE]
  }
  if (nrow(rows) == 0) {
    return("  (none)")
  }
  rows <- head(rows, n)
  paste0(
    "  ",
    seq_len(nrow(rows)), ". ",
    rows$phenotype, " ~ ", rows$predictor,
    ": ", col, " = ", sprintf("%.3f", rows[[col]]),
    ", p = ", sprintf("%.3g", rows[[p_col]]),
    collapse = "\n"
  )
}

build_summary <- function(cors, dat) {
  sign_flips <- cors |>
    filter(!is.na(pearson_r), !is.na(spearman_rho)) |>
    filter(sign(pearson_r) != sign(spearman_rho)) |>
    filter(abs(pearson_r) > SIGN_FLIP_THRESHOLD | abs(spearman_rho) > SIGN_FLIP_THRESHOLD) |>
    arrange(desc(abs(pearson_r - spearman_rho)))

  nonlinear <- bind_rows(lapply(PHENOTYPES, function(ph) {
    bind_rows(lapply(CONTINUOUS_PREDICTORS, function(pr) {
      data.frame(
        phenotype = ph,
        predictor = pr,
        nonlinear = detect_nonlinearity(dat, ph, pr),
        stringsAsFactors = FALSE
      )
    }))
  })) |>
    filter(nonlinear)

  outlier_rows <- bind_rows(lapply(PHENOTYPES, function(ph) {
    bind_rows(lapply(PREDICTORS, function(pr) {
      sp <- detect_outliers(dat, ph, pr)
      if (length(sp) == 0) {
        return(NULL)
      }
      data.frame(
        phenotype = ph,
        predictor = pr,
        species = sp,
        stringsAsFactors = FALSE
      )
    }))
  }))

  large_gap <- cors |>
    filter(!is.na(pearson_r), !is.na(spearman_rho)) |>
    mutate(gap = abs(pearson_r - spearman_rho)) |>
    filter(gap >= NONLINEARITY_THRESHOLD) |>
    arrange(desc(gap))

  lines <- c(
    "Exploratory trait–innovation summary",
    paste0("Species in table: ", nrow(dat)),
    paste0("Phenotypes: ", paste(PHENOTYPES, collapse = ", ")),
    "",
    "Top positive Pearson correlations",
    format_top_list(cors, "pearson_r", "positive"),
    "",
    "Top negative Pearson correlations",
    format_top_list(cors, "pearson_r", "negative"),
    "",
    "Top positive Spearman correlations",
    format_top_list(cors, "spearman_rho", "positive", p_col = "spearman_p"),
    "",
    "Top negative Spearman correlations",
    format_top_list(cors, "spearman_rho", "negative", p_col = "spearman_p"),
    "",
    "Predictors with sign flip between Pearson and Spearman",
    if (nrow(sign_flips) == 0) {
      "  (none above threshold)"
    } else {
      paste0(
        "  ",
        sign_flips$phenotype, " ~ ", sign_flips$predictor,
        ": r = ", sprintf("%.3f", sign_flips$pearson_r),
        ", rho = ", sprintf("%.3f", sign_flips$spearman_rho),
        collapse = "\n"
      )
    },
    "",
    "Predictors with strong nonlinearity (quadratic improvement)",
    if (nrow(nonlinear) == 0) {
      "  (none detected)"
    } else {
      paste0("  ", nonlinear$phenotype, " ~ ", nonlinear$predictor, collapse = "\n")
    },
    "",
    "Large Pearson vs Spearman gaps (>= 0.15)",
    if (nrow(large_gap) == 0) {
      "  (none)"
    } else {
      paste0(
        "  ",
        head(large_gap, 10)$phenotype, " ~ ", head(large_gap, 10)$predictor,
        ": gap = ", sprintf("%.3f", head(large_gap, 10)$gap),
        collapse = "\n"
      )
    },
    "",
    "Possible outliers (|standardized residual| > 2.5)",
    if (nrow(outlier_rows) == 0) {
      "  (none flagged)"
    } else {
      paste0(
        "  ",
        outlier_rows$phenotype, " ~ ", outlier_rows$predictor, ": ",
        outlier_rows$species,
        collapse = "\n"
      )
    }
  )
  paste(lines, collapse = "\n")
}

# ---- main ----
main <- function() {
  dir_ok(OUT_DIR)
  dir_ok(FIG_DIR)

  dat <- load_cleaned_data()

  cors <- compute_all_correlations(dat)
  write_tsv(cors, file.path(OUT_DIR, "correlations.tsv"))
  message("Wrote correlations.tsv (", nrow(cors), " rows)")

  for (ph in PHENOTYPES) {
    ph_fig_dir <- file.path(FIG_DIR, ph)
    dir_ok(ph_fig_dir)

    for (pr in CONTINUOUS_PREDICTORS) {
      stats_row <- cors[cors$phenotype == ph & cors$predictor == pr, , drop = FALSE]
      plot_continuous_scatter(
        dat, ph, pr, stats_row,
        file.path(ph_fig_dir, paste0("scatter_", pr, ".pdf"))
      )
    }

    plot_dietbreadth(dat, ph, ph_fig_dir)

    for (pr in CATEGORICAL_PREDICTORS) {
      plot_categorical(dat, ph, pr, ph_fig_dir)
    }
  }
  message("Wrote scatter and group plots to ", FIG_DIR)

  plot_correlation_matrix(
    dat, CORRELATION_MATRIX_TRAITS,
    file.path(FIG_DIR, "correlation_matrix.pdf")
  )
  message("Wrote correlation matrix")

  summary_text <- build_summary(cors, dat)
  writeLines(summary_text, file.path(OUT_DIR, "summary.txt"))
  cat("\n", summary_text, "\n", sep = "")

  message("Done. Output in: ", OUT_DIR)
}

if (sys.nframe() == 0) {
  main()
}
