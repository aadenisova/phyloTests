#!/usr/bin/env Rscript
# contMap-style semicircle plot for TOTALINNOVATIONS2025_ResEff with PhyloPic silhouettes.

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(rphylopic)
})

PHENOTYPE <- "TOTALINNOVATIONS2025_ResEff"
ROOT <- normalizePath(".", wins = FALSE)
BASE <- file.path(ROOT, "results", PHENOTYPE)
# OUT <- file.path(BASE, "figures", "contMap_fan_phylopic.pdf")
OUT <- file.path(BASE, "figures", "contMap_fan_phylopic.png")

CACHE <- file.path(ROOT, "data", "phylopic_cache.rds")
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

get_uuid_safe <- function(name) {
  for (n in c(name, strsplit(name, " ", fixed = TRUE)[[1]][1])) {
    u <- tryCatch(get_uuid(n, n = 1), error = function(e) NA_character_)
    if (!is.na(u)) return(u)
  }
  NA_character_
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

cm <- contMap(tree, x, plot = FALSE, lut = c("blue", "white", "red"))
# pdf(OUT, width = 14, height = 8)
png(
  OUT,
  width = 4000,
  height = 2400,
  res = 300,
  bg = "transparent"
)
par(mar = c(1, 1, 2, 1), xpd = NA)
plot(cm, type = "arc", fsize = 0.7, outline = FALSE, lwd = 3)

pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
img_h <- 0.015

img_scale <- c(
  "Sarcoramphus papa"   = 0.65,
  "Gypaetus barbatus"   = 0.65,
  "Aquila chrysaetos"   = 0.65,
  "Anas platyrhynchos"  = 0.70,
  "Clangula hyemalis"   = 0.70,
  "Bucephala clangula"  = 0.70,
  "Gallinula chloropus" = 0.70
)

for (t in seq_len(Ntip(tree))) {
  if (is.null(imgs[[t]])) next
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
  # pos <- xy + dir * (strwidth(tt, cex = pp$cex) + img_h * 0.8)
  h <- img_h
  if (lab %in% names(img_scale))
      h <- img_h * img_scale[lab]

  pos <- xy + dir * (strwidth(tt, cex = pp$cex) + h * 0.8)

  # add_phylopic_base(img = imgs[[t]], x = pos[1], y = pos[2], height = img_h)
  

  if (lab %in% names(img_scale))
      h <- img_h * img_scale[lab]

  add_phylopic_base(
      img = imgs[[t]],
      x = pos[1],
      y = pos[2],
      height = h
  )
}

dev.off()
message("Saved: ", OUT)
