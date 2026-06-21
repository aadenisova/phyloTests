#!/usr/bin/env Rscript
# contMap semicircle with nested foreground-branch arc tree and PhyloPic silhouettes.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(rphylopic)
})

PHENOTYPE <- "TOTALINNOVATIONS2025_ResEff"
ROOT <- normalizePath(".", wins = FALSE)
BASE <- file.path(ROOT, "results", PHENOTYPE)
OUT <- file.path(BASE, "figures", "contMap_foreground_fan_phylopic.pdf")
CACHE <- file.path(ROOT, "data", "phylopic_cache.rds")
INNER_SCALE <- 0.35
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

get_uuid_safe <- function(name) {
  for (n in c(name, strsplit(name, " ", fixed = TRUE)[[1]][1])) {
    u <- tryCatch(get_uuid(n, n = 1), error = function(e) NA_character_)
    if (!is.na(u)) return(u)
  }
  NA_character_
}

tip_img_pos <- function(pp, t, img_h, label_gap = 0.012) {
  xy <- c(pp$xx[as.character(t)], pp$yy[as.character(t)])
  aa <- atan2(xy[2], xy[1])
  dir <- c(cos(aa), sin(aa))
  xy + dir * (label_gap + img_h * 0.7)
}

tree <- read.tree(file.path(ROOT, "data", "roadies_birds_allbirdtraits.nwk"))
xdf <- read.delim(file.path(BASE, "continuous_phenotype.tsv"), stringsAsFactors = FALSE)
tree <- keep.tip(tree, xdf$species)
x <- setNames(xdf$innovation, xdf$species)[tree$tip.label]
binomial <- gsub("_", " ", tree$tip.label)

cache <- if (file.exists(CACHE)) readRDS(CACHE) else list()
for (nm in unique(binomial)) {
  if (!is.null(cache[[nm]])) next
  uuid <- get_uuid_safe(nm)
  if (!is.na(uuid)) cache[[nm]] <- get_phylopic(uuid = uuid)
}
saveRDS(cache, CACHE)
imgs <- lapply(binomial, function(nm) cache[[nm]])

foreground <- read.delim(file.path(BASE, "foreground_branches.tsv"), stringsAsFactors = FALSE)
branch_changes <- read.delim(file.path(BASE, "branch_changes.tsv"), stringsAsFactors = FALSE)
ne <- Nedge(tree)
edge_col <- rep("gray60", ne)
edge_w <- rep(0.6, ne)
idx <- match(foreground$branch_id, branch_changes$branch_id)
idx <- idx[!is.na(idx)]
edge_col[idx] <- "red"
edge_w[idx] <- 2.5

inner <- tree
inner$edge.length <- inner$edge.length * INNER_SCALE

cm <- contMap(tree, x, plot = FALSE, lut = c("blue", "white", "red"))
pdf(OUT, width = 14, height = 8)
par(mar = c(1, 1, 2, 1), xpd = NA)
plot(cm, type = "arc", fsize = 0.7, outline = FALSE, lwd = 3, ftype = "off")
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
plotTree(
  inner, type = "arc", add = TRUE,
  edge.color = edge_col, edge.width = edge_w,
  ftype = "off", fsize = 0.5
)
title("Total innovation (research-adjusted)", cex.main = 1.2)
img_h <- 0.02
for (t in seq_len(Ntip(tree))) {
  if (is.null(imgs[[t]])) next
  pos <- tip_img_pos(pp, t, img_h)
  add_phylopic_base(img = imgs[[t]], x = pos[1], y = pos[2], height = img_h)
}
dev.off()
message("Saved: ", OUT)
