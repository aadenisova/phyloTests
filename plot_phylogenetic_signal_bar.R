#!/usr/bin/env Rscript
# Bar plot of Blomberg's K from phylogenetic_signal_master.tsv.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

PROJECT_ROOT <- normalizePath(".", wins = FALSE)
MASTER_FILE <- file.path(PROJECT_ROOT, "results", "_master", "phylogenetic_signal_master.tsv")
OUT_DIR <- file.path(PROJECT_ROOT, "final_figure")

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

BASE_SIZE <- 14
TITLE_SIZE <- 18
ALPHA <- 0.05

labelize <- function(x) {
  out <- unname(TRAIT_LABELS[x])
  ifelse(is.na(out), gsub("_", " ", x), out)
}

df <- read.delim(MASTER_FILE, stringsAsFactors = FALSE)
df$K_p_adj <- p.adjust(df$K_p, method = "BH")
df$trait_label <- factor(
  labelize(df$trait),
  levels = labelize(df$trait[order(df$K, decreasing = TRUE)])
)
df$significant <- df$K_p_adj <= ALPHA

p <- ggplot(df, aes(x = trait_label, y = K)) +
  geom_col(fill = "#4575B4", width = 0.72) +
  geom_text(
    data = df[df$significant, , drop = FALSE],
    aes(label = "*"),
    vjust = -0.4,
    size = 6,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Phylogenetic signal (Blomberg's K)",
    x = NULL,
    y = expression(italic(K))
  ) +
  theme_bw(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(size = TITLE_SIZE, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE),
    axis.text.y = element_text(size = BASE_SIZE),
    axis.title.y = element_text(size = BASE_SIZE + 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(OUT_DIR, "phylogenetic_signal_bar.pdf"), p, width = 10, height = 6)
ggsave(file.path(OUT_DIR, "phylogenetic_signal_bar.png"), p, width = 10, height = 6, dpi = 300)

message("Saved: ", file.path(OUT_DIR, "phylogenetic_signal_bar.pdf"))
message("Saved: ", file.path(OUT_DIR, "phylogenetic_signal_bar.png"))
