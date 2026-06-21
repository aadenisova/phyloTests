#!/usr/bin/env Rscript
# Foreground-branch semicircle tree + innovation ring + PhyloPic silhouettes.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(rphylopic)
})

PHENOTYPE <- "TOTALINNOVATIONS2025_ResEff"
ROOT <- normalizePath(".", wins = FALSE)
BASE <- file.path(ROOT, "results", PHENOTYPE)
OUT <- file.path(BASE, "figures", "foreground_fan_phylopic.pdf")
CACHE <- file.path(ROOT, "data", "phylopic_cache.rds")
LUT <- c("blue", "white", "red")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

tip_colors <- function(values, lut = LUT) {
  norm <- (values - min(values)) / diff(range(values))
  rgb(colorRamp(lut)(norm), maxColorValue = 255)
}

draw_inno_ring <- function(pp, tip_labels, values, r1, r2) {
  n <- length(tip_labels)
  ang <- atan2(pp$yy[as.character(seq_len(n))], pp$xx[as.character(seq_len(n))])
  ord <- order(ang)
  ang <- ang[ord]
  cols <- tip_colors(values[tip_labels[ord]])
  d <- diff(ang)
  br <- c(ang[1] - d[1] / 2, (head(ang, -1) + tail(ang, -1)) / 2, ang[n] + tail(d, 1) / 2)
  for (i in seq_len(n)) {
    th <- seq(br[i], br[i + 1], length.out = 40)
    polygon(
      c(r1 * cos(th), rev(r2 * cos(th))),
      c(r1 * sin(th), rev(r2 * sin(th))),
      col = cols[i], border = NA
    )
  }
}

tip_img_pos <- function(pp, t, img_h, ring_w) {
  xy <- c(pp$xx[as.character(t)], pp$yy[as.character(t)])
  lab <- binomial[t]
  aa <- atan2(xy[2], xy[1]) * 180 / pi
  tt <- if (aa > 90 && aa < 270) {
    paste(lab, paste(rep(" ", pp$label.offset), collapse = ""), sep = "")
  } else {
    paste(paste(rep(" ", pp$label.offset), collapse = ""), lab, sep = "")
  }
  srt <- if (aa > 90 && aa < 270) 180 + aa else aa
  dir <- c(cos(srt * pi / 180), sin(srt * pi / 180))
  if (aa > 90 && aa < 270) dir <- -dir
  xy + dir * (strwidth(tt, cex = pp$cex) + ring_w + img_h * 0.8)
}

tree <- read.tree(file.path(ROOT, "data", "roadies_birds_allbirdtraits.nwk"))
xdf <- read.delim(file.path(BASE, "continuous_phenotype.tsv"), stringsAsFactors = FALSE)
tree <- keep.tip(tree, xdf$species)
x <- setNames(xdf$innovation, xdf$species)[tree$tip.label]
binomial <- gsub("_", " ", tree$tip.label)
cache <- if (file.exists(CACHE)) readRDS(CACHE) else list()
imgs <- lapply(binomial, function(nm) cache[[nm]])

ne <- Nedge(tree)
edge_col <- rep("gray60", ne)
edge_w <- rep(0.6, ne)
foreground <- read.delim(file.path(BASE, "foreground_branches.tsv"), stringsAsFactors = FALSE)
branch_changes <- read.delim(file.path(BASE, "branch_changes.tsv"), stringsAsFactors = FALSE)
idx <- match(foreground$branch_id, branch_changes$branch_id)
idx <- idx[!is.na(idx)]
edge_col[idx] <- "red"
edge_w[idx] <- 2.5

pdf(OUT, width = 14, height = 8)
par(mar = c(1, 1, 2, 1), xpd = NA)
plotTree(tree, type = "arc", edge.color = edge_col, edge.width = edge_w, ftype = "off", fsize = 0.7)
plotTree(
  tree, type = "arc", edge.color = rep("transparent", ne), edge.width = edge_w,
  ftype = "i", fsize = 0.7, add = TRUE
)
title("Foreground branches (independent innovation candidates)", cex.main = 1.2)

pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
label_r <- max(vapply(seq_len(Ntip(tree)), function(t) {
  xy <- c(pp$xx[as.character(t)], pp$yy[as.character(t)])
  sqrt(sum(xy^2)) + strwidth(binomial[t], cex = pp$cex)
}, numeric(1)))
r1 <- label_r * 1.03
r2 <- label_r * 1.10
draw_inno_ring(pp, tree$tip.label, x, r1, r2)

img_h <- 0.02
ring_w <- r2 - r1
for (t in seq_len(Ntip(tree))) {
  if (is.null(imgs[[t]])) next
  pos <- tip_img_pos(pp, t, img_h, ring_w)
  add_phylopic_base(img = imgs[[t]], x = pos[1], y = pos[2], height = img_h)
}
dev.off()
message("Saved: ", OUT)
