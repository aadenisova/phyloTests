#!/usr/bin/env Rscript
# Build publication multi-panel figure and formatted Excel workbook.

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(openxlsx)
})

PROJECT_ROOT <- normalizePath(".", wins = FALSE)
RESULTS_DIR <- file.path(PROJECT_ROOT, "results")
MASTER_DIR <- file.path(RESULTS_DIR, "_master")
PAPER_DIR <- file.path(RESULTS_DIR, "paper")
TREE_FILE <- file.path(PROJECT_ROOT, "data", "tree_clootl.newick")
IMAGE_DIR <- file.path(PROJECT_ROOT, "data", "species_images")

PHENOTYPE_MAP <- c(
  TOTAL = "TOTALINNOVATIONS2025_ResEff",
  FOOD = "FOODINNO2025_ResEff",
  TECH = "TECHINNO2025_ResEff"
)
PHENOTYPE_LABELS <- c(
  TOTAL = "Total innovation",
  FOOD = "Food innovation",
  TECH = "Technical innovation"
)
INNOVATION_TRAITS <- unname(PHENOTYPE_MAP)
ALPHA <- 0.05
BASE_SIZE <- 14
TITLE_SIZE <- 18
SUBTITLE_SIZE <- 13
TAG_SIZE <- 16
SIG_COLOR <- "#2CA25F"
NS_COLOR <- "#888888"

TRAIT_LABELS <- c(
  Brain_size = "Brain size",
  Relative_brain_size = "Relative brain size",
  Mass = "Body mass",
  GenerationLength = "Generation length",
  Range.Size = "Range size",
  HabitatBreadth = "Habitat breadth",
  DietBreadth = "Diet breadth",
  UrbanFULL = "Urban tolerance",
  Migration = "Migration",
  Trophic_level = "Trophic level",
  TOTALINNOVATIONS2025_ResEff = "Total innovation",
  FOODINNO2025_ResEff = "Food innovation",
  TECHINNO2025_ResEff = "Technical innovation"
)

dir_ok <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)

labelize <- function(x) {
  out <- unname(TRAIT_LABELS[x])
  ifelse(is.na(out), gsub("_", " ", x), out)
}

labelize_phenotype <- function(x) {
  out <- unname(PHENOTYPE_LABELS[x])
  ifelse(is.na(out), x, out)
}

value_to_color <- function(values, lut = c("#4575B4", "#FFFFFF", "#D73027")) {
  vals <- values[!is.na(values)]
  rng <- range(vals)
  if (length(rng) < 2 || diff(rng) == 0) rng[2] <- rng[1] + 1e-9
  norm <- (values - rng[1]) / diff(rng)
  norm[is.na(norm)] <- 0.5
  rgb(grDevices::colorRamp(lut)(norm), maxColorValue = 255)
}

load_analysis_tree <- function() {
  species <- Reduce(
    intersect,
    lapply(PHENOTYPE_MAP, function(folder) {
      read.delim(
        file.path(RESULTS_DIR, folder, "continuous_phenotype.tsv"),
        stringsAsFactors = FALSE
      )$species
    })
  )
  tree <- read.tree(TREE_FILE)
  keep.tip(tree, intersect(species, tree$tip.label))
}

load_phenotype_data <- function(code) {
  folder <- PHENOTYPE_MAP[[code]]
  base <- file.path(RESULTS_DIR, folder)
  list(
    code = code,
    label = PHENOTYPE_LABELS[[code]],
    folder = folder,
    x = setNames(
      read.delim(file.path(base, "continuous_phenotype.tsv"), stringsAsFactors = FALSE)$innovation,
      read.delim(file.path(base, "continuous_phenotype.tsv"), stringsAsFactors = FALSE)$species
    ),
    branch_changes = read.delim(file.path(base, "branch_changes.tsv"), stringsAsFactors = FALSE),
    foreground = read.delim(file.path(base, "foreground_branches.tsv"), stringsAsFactors = FALSE)
  )
}

plot_tree_png <- function(tree, branch_changes, foreground, x, title, out_png) {
  ne <- Nedge(tree)
  edge_col <- rep("gray70", ne)
  edge_w <- rep(0.8, ne)
  idx <- match(foreground$branch_id, branch_changes$branch_id)
  idx <- idx[!is.na(idx)]
  edge_col[idx] <- "#E41A1C"
  edge_w[idx] <- 2.2
  tip_col <- value_to_color(as.numeric(x[tree$tip.label]))

  png(out_png, width = 1400, height = 1400, res = 160, bg = "white")
  par(mar = c(0.5, 0.5, 2.2, 0.5), xpd = NA)
  plot(
    tree,
    type = "fan",
    show.tip.label = FALSE,
    edge.color = edge_col,
    edge.width = edge_w,
    tip.color = tip_col,
    no.margin = TRUE
  )
  title(main = title, cex.main = 1.35, font.main = 2)
  legend(
    "bottomleft",
    legend = c("Foreground branch", "Tip innovativeness"),
    col = c("#E41A1C", "#4575B4"),
    lwd = c(2.5, 5),
    bty = "n",
    cex = 1.0
  )
  dev.off()
}

make_tree_panel <- function(tree, pdata) {
  tmp <- tempfile(fileext = ".png")
  x <- pdata$x[tree$tip.label]
  plot_tree_png(tree, pdata$branch_changes, pdata$foreground, x, pdata$label, tmp)
  img <- png::readPNG(tmp)
  ggdraw() +
    draw_grob(
      grid::rasterGrob(img, width = grid::unit(1, "npc"), height = grid::unit(1, "npc"), interpolate = TRUE)
    )
}

make_forest_plot <- function(df) {
  df <- df |>
    mutate(
      trait_label = labelize(trait),
      phenotype_label = factor(
        labelize_phenotype(phenotype),
        levels = labelize_phenotype(names(PHENOTYPE_MAP))
      ),
      significant = p < ALPHA,
      lower = beta - 1.96 * se,
      upper = beta + 1.96 * se
    )

  ggplot(df, aes(x = beta, y = trait_label, color = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.3, linewidth = 0.7) +
    geom_point(size = 3.2) +
    facet_wrap(~phenotype_label, nrow = 1) +
    scale_color_manual(
      values = c("FALSE" = NS_COLOR, "TRUE" = SIG_COLOR),
      labels = c("FALSE" = expression(p >= 0.05), "TRUE" = expression(p < 0.05))
    ) +
    labs(
      title = "Univariate PGLS: ecological predictors",
      x = expression(beta ~ "(z-scaled predictors)"),
      y = NULL,
      color = NULL
    ) +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(size = TITLE_SIZE, face = "bold"),
      axis.text = element_text(size = BASE_SIZE),
      axis.title = element_text(size = BASE_SIZE + 1),
      legend.text = element_text(size = BASE_SIZE),
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold", size = BASE_SIZE)
    )
}

make_parent_child_plot <- function(codes, tree) {
  rows <- lapply(codes, function(code) {
    pd <- load_phenotype_data(code)
    bc <- pd$branch_changes
    bc$foreground <- bc$branch_id %in% pd$foreground$branch_id
    bc$phenotype <- code
    bc$phenotype_label <- factor(
      PHENOTYPE_LABELS[[code]],
      levels = labelize_phenotype(names(PHENOTYPE_MAP))
    )
    bc$point_label <- ifelse(
      bc$foreground,
      sprintf("\u0394=%.2f", bc$delta),
      NA_character_
    )
    bc
  })
  df <- bind_rows(rows)
  fg <- df[df$foreground, , drop = FALSE]

  ggplot(df, aes(parent_value, child_value)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.6) +
    geom_point(
      data = df[!df$foreground, , drop = FALSE],
      aes(color = foreground),
      shape = 16, alpha = 0.35, size = 2.2
    ) +
    geom_point(
      data = fg,
      aes(color = foreground),
      shape = 17, alpha = 0.95, size = 3.8
    ) +
    geom_text(
      data = fg,
      aes(label = point_label),
      size = BASE_SIZE / 3.2,
      fontface = "bold",
      color = "#B2182B",
      nudge_y = 0.08,
      check_overlap = TRUE
    ) +
    facet_wrap(~phenotype_label, nrow = 1) +
    scale_color_manual(
      values = c("FALSE" = NS_COLOR, "TRUE" = "#E41A1C"),
      labels = c("Background branch", "Foreground branch")
    ) +
    labs(
      title = "Parent vs child values along branches",
      x = "Parent value (ancestral / tip)",
      y = "Child value (descendant / tip)",
      color = NULL
    ) +
    theme_bw(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(size = TITLE_SIZE, face = "bold"),
      axis.text = element_text(size = BASE_SIZE),
      axis.title = element_text(size = BASE_SIZE + 1),
      legend.text = element_text(size = BASE_SIZE),
      legend.position = "bottom",
      strip.background = element_rect(fill = "gray95"),
      strip.text = element_text(face = "bold", size = BASE_SIZE)
    )
}

build_figure <- function(tree) {
  pdata <- lapply(names(PHENOTYPE_MAP), load_phenotype_data)
  names(pdata) <- names(PHENOTYPE_MAP)

  p_trees <- lapply(names(PHENOTYPE_MAP), function(code) {
    make_tree_panel(tree, pdata[[code]])
  })

  pgls_uni <- read.delim(
    file.path(MASTER_DIR, "pgls_univariate_master.tsv"),
    stringsAsFactors = FALSE
  )
  p_forest <- make_forest_plot(pgls_uni)
  p_parent <- make_parent_child_plot(names(PHENOTYPE_MAP), tree)

  panel_trees <- wrap_plots(
    p_trees[[1]] + labs(tag = "A"),
    p_trees[[2]] + labs(tag = "B"),
    p_trees[[3]] + labs(tag = "C"),
    nrow = 1
  )

  fig <- panel_trees / (p_forest + labs(tag = "D")) / (p_parent + labs(tag = "E")) +
    plot_layout(heights = c(2, 0.55, 0.45)) +
    plot_annotation(
      title = "Phylogenetic patterns of bird innovativeness",
      subtitle = paste0(
        "Research-adjusted innovation scores (n = ", Ntip(tree),
        " species). Red branches = independent innovation candidates (\u0394 > mean + 1 SD)."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = TITLE_SIZE + 2),
        plot.subtitle = element_text(size = SUBTITLE_SIZE, color = "gray30"),
        plot.tag = element_text(size = TAG_SIZE, face = "bold")
      )
    )

  fig
}

humanize_master <- function(df, type = c("signal", "uni", "mv", "top20")) {
  type <- match.arg(type)
  switch(type,
    signal = df |>
      mutate(
        Trait = labelize(trait),
        Phenotype = labelize_phenotype(phenotype)
      ) |>
      select(Trait, Phenotype, K, `K p-value` = K_p, Lambda = lambda, `Lambda p-value` = lambda_p),
    uni = df |>
      mutate(
        Predictor = labelize(trait),
        Phenotype = labelize_phenotype(phenotype),
        `Adj. R2` = r2
      ) |>
      select(Phenotype, Predictor, beta, se, p, `Adj. R2`),
    mv = df |>
      mutate(
        Predictor = labelize(trait),
        Phenotype = labelize_phenotype(phenotype),
        `Adj. R2` = r2
      ) |>
      select(Phenotype, Predictor, beta, se, p, `Adj. R2`),
    top20 = df |>
      mutate(Phenotype = labelize_phenotype(phenotype)) |>
      select(
        Phenotype, `Branch ID` = branch_id, Parent = parent, Child = child,
        `Parent value` = parent_value, `Child value` = child_value,
        `Branch length` = branch_length, Delta = delta
      )
  )
}

style_sheet <- function(wb, sheet, df, header_color = "#1F4E79") {
  header_style <- createStyle(
    fontColour = "#FFFFFF",
    fgFill = header_color,
    halign = "center",
    textDecoration = "bold",
    border = "Bottom",
    borderColour = "#CCCCCC"
  )
  body_style <- createStyle(border = "Bottom", borderColour = "#E6E6E6")
  addStyle(wb, sheet, header_style, rows = 1, cols = seq_len(ncol(df)), gridExpand = TRUE)
  if (nrow(df) > 0) {
    addStyle(wb, sheet, body_style, rows = 2:(nrow(df) + 1), cols = seq_len(ncol(df)), gridExpand = TRUE)
  }
  setColWidths(wb, sheet, cols = seq_len(ncol(df)), widths = "auto")
  freezePane(wb, sheet, firstActiveRow = 2)
}

build_species_sheet <- function(tree) {
  wide <- Reduce(
    function(x, y) merge(x, y, by = "species", all = TRUE),
    lapply(names(PHENOTYPE_MAP), function(code) {
      df <- read.delim(
        file.path(RESULTS_DIR, PHENOTYPE_MAP[[code]], "continuous_phenotype.tsv"),
        stringsAsFactors = FALSE
      )
      names(df)[2] <- code
      df
    })
  ) |>
    filter(species %in% tree$tip.label) |>
    mutate(`Scientific name` = gsub("_", " ", species)) |>
    select(
      `Scientific name`, Species = species,
      `Total innovation` = TOTAL,
      `Food innovation` = FOOD,
      `Technical innovation` = TECH
    ) |>
    arrange(`Scientific name`)

  if (dir.exists(IMAGE_DIR)) {
    wide$Image <- vapply(wide$Species, function(sp) {
      hits <- list.files(IMAGE_DIR, pattern = paste0("^", sp, "\\.(jpg|jpeg|png|webp)$"), ignore.case = TRUE)
      if (length(hits)) hits[1] else NA_character_
    }, character(1))
  }

  wide
}

insert_species_images <- function(wb, sheet, df) {
  if (!"Image" %in% names(df)) return(invisible(NULL))
  img_col <- which(names(df) == "Image")
  for (i in seq_len(nrow(df))) {
    img_file <- df$Image[i]
    if (is.na(img_file) || !nzchar(img_file)) next
    path <- file.path(IMAGE_DIR, img_file)
    if (!file.exists(path)) next
    tryCatch(
      insertImage(
        wb, sheet, path,
        startRow = i + 1, startCol = img_col,
        width = 1.8, height = 1.8, units = "cm"
      ),
      error = function(e) NULL
    )
  }
}

build_excel <- function(tree) {
  signal <- read.delim(file.path(MASTER_DIR, "phylogenetic_signal_master.tsv"), stringsAsFactors = FALSE)
  uni <- read.delim(file.path(MASTER_DIR, "pgls_univariate_master.tsv"), stringsAsFactors = FALSE)
  mv <- read.delim(file.path(MASTER_DIR, "pgls_multivariate_master.tsv"), stringsAsFactors = FALSE)
  top20 <- read.delim(file.path(MASTER_DIR, "top20_branch_changes_master.tsv"), stringsAsFactors = FALSE)
  beta_mat <- read.delim(file.path(MASTER_DIR, "pgls_univariate_beta_matrix.tsv"), stringsAsFactors = FALSE)
  lambda_mat <- read.delim(file.path(MASTER_DIR, "innovation_lambda_matrix.tsv"), stringsAsFactors = FALSE)

  overview <- data.frame(
    Item = c(
      "Analysis", "Tree source", "Species included",
      "Innovation phenotypes", "Foreground threshold", "Significance level"
    ),
    Value = c(
      "Phylogenetic comparative analysis of bird innovativeness",
      "Clootl bird tree",
      as.character(Ntip(tree)),
      "Total, Food, Technical (research-adjusted)",
      "Branch delta > mean + 1 SD",
      paste0("p < ", ALPHA)
    ),
    stringsAsFactors = FALSE
  )

  species_df <- build_species_sheet(tree)
  beta_mat$Predictor <- labelize(beta_mat$trait)
  beta_mat$trait <- NULL
  beta_mat <- beta_mat |> select(Predictor, everything())

  wb <- createWorkbook()
  addWorksheet(wb, "Overview")
  addWorksheet(wb, "Species")
  addWorksheet(wb, "Phylogenetic signal")
  addWorksheet(wb, "PGLS univariate")
  addWorksheet(wb, "PGLS multivariate")
  addWorksheet(wb, "Top branch changes")
  addWorksheet(wb, "Beta matrix")
  addWorksheet(wb, "Lambda matrix")

  writeData(wb, "Overview", overview)
  writeData(wb, "Species", species_df)
  writeData(wb, "Phylogenetic signal", humanize_master(signal, "signal"))
  writeData(wb, "PGLS univariate", humanize_master(uni, "uni"))
  writeData(wb, "PGLS multivariate", humanize_master(mv, "mv"))
  writeData(wb, "Top branch changes", humanize_master(top20, "top20"))
  writeData(wb, "Beta matrix", beta_mat)
  writeData(wb, "Lambda matrix", lambda_mat)

  for (sh in names(wb$worksheets)) {
    df <- switch(sh,
      Overview = overview,
      Species = species_df,
      `Phylogenetic signal` = humanize_master(signal, "signal"),
      `PGLS univariate` = humanize_master(uni, "uni"),
      `PGLS multivariate` = humanize_master(mv, "mv"),
      `Top branch changes` = humanize_master(top20, "top20"),
      `Beta matrix` = beta_mat,
      `Lambda matrix` = lambda_mat
    )
    style_sheet(wb, sh, df)
  }

  insert_species_images(wb, "Species", species_df)
  saveWorkbook(wb, file.path(PAPER_DIR, "phylogenetic_innovation_results.xlsx"), overwrite = TRUE)
}

main <- function() {
  dir_ok(PAPER_DIR)
  tree <- load_analysis_tree()
  message("Building figure for ", Ntip(tree), " species...")

  fig <- build_figure(tree)
  ggsave(
    file.path(PAPER_DIR, "phylogenetic_innovation_figure.pdf"),
    fig, width = 18, height = 20, limitsize = FALSE
  )
  ggsave(
    file.path(PAPER_DIR, "phylogenetic_innovation_figure.png"),
    fig, width = 18, height = 20, dpi = 300, limitsize = FALSE
  )

  message("Building Excel workbook...")
  build_excel(tree)

  message("Done:")
  message("  ", file.path(PAPER_DIR, "phylogenetic_innovation_figure.pdf"))
  message("  ", file.path(PAPER_DIR, "phylogenetic_innovation_figure.png"))
  message("  ", file.path(PAPER_DIR, "phylogenetic_innovation_results.xlsx"))
}

if (sys.nframe() == 0) {
  main()
}
