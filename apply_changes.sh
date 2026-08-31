#!/bin/bash
set -e
cd "$(dirname "$0")"

git rm -q --cached config.example.yml 2>/dev/null || true
rm -f config.example.yml config.yml

cat > .gitignore <<'GEOF'
# inputs: download from GEO (GSE295892, GSE344451) — see README
P560-2/data/
P560-2/DoRothEA/
P560-4/data/
P560-4/DoRothEA/

# generated
P560-2/results/
P560-4/results/
shared_results/
docs/

*.RData
*.rds
.Rhistory
.RData
.Rproj.user/
.DS_Store
GEOF

cat > R/00_setup.R <<'GEOF'
# Sourced by every notebook. Set DATASET first.

library(edgeR); library(limma); library(tidyverse); library(ggrepel)
library(org.Hs.eg.db); library(clusterProfiler); library(enrichplot); library(cowplot)
library(svglite)

library(conflicted)
for (fn in c("select","filter","rename","slice","mutate","count","desc","first"))
  conflict_prefer(fn, "dplyr", quiet = TRUE)
for (fn in c("setdiff","intersect","union","setequal"))
  conflict_prefer(fn, "base", quiet = TRUE)
conflict_prefer("simplify", "clusterProfiler", quiet = TRUE)

DATASETS <- list(
  "P560-2" = list(
    label       = "Dynasore",
    inhibitor   = "Dynasore",
    vehicle     = "DMSO",
    base_arm    = "Base_Dynasore",
    counts      = "240201-P560-2_RawGeneCounts.csv",
    interaction = TRUE
  ),
  "P560-4" = list(
    label       = "panJAKi",
    inhibitor   = "JAKi",
    vehicle     = "DMSO",
    base_arm    = NA_character_,   # no JAKi-only arm; interaction not identifiable
    counts      = "P560-4_gene_counts.csv",
    interaction = FALSE
  )
)

if (!exists("DATASET")) stop("set DATASET before sourcing 00_setup.R")
if (!DATASET %in% names(DATASETS)) stop("unknown DATASET: ", DATASET)

CFG <- DATASETS[[DATASET]]

# inputs live under <study code>/, as documented in the README
STUDY_DIR   <- DATASET
COUNT_FILE  <- file.path(STUDY_DIR, "data", CFG$counts)
REGULON_FILE <- file.path(STUDY_DIR, "DoRothEA", "BEST_viperRegulon.rdata")

missing <- c(COUNT_FILE, REGULON_FILE)[!file.exists(c(COUNT_FILE, REGULON_FILE))]
if (length(missing))
  stop("missing input files:\n  ", paste(missing, collapse = "\n  "),
       "\nsee README section 4")

OUT_DIR    <- file.path(STUDY_DIR, "results")
FIG_DIR    <- file.path(OUT_DIR, "figures")
DATA_DIR   <- file.path(OUT_DIR, "tables")
STATE_DIR  <- file.path(OUT_DIR, "state")
SHARED_DIR <- file.path("shared_results")
for (d in c(OUT_DIR, FIG_DIR, DATA_DIR, STATE_DIR, SHARED_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

STATE_FILE <- file.path(STATE_DIR, "analysis_state.RData")

LFC_CUTOFF <- 1      # log2 scale
FDR_CUTOFF <- 0.05
TREAT_LFC  <- 0
N_LABEL    <- 20
ANNO_COLS  <- 4      # geneId, geneName, description, geneBiotype

CYTOKINES  <- c("IL19","IL20","IL22","IL24","IL26","IFNL3")
CYT_LABEL  <- c(IL19="IL-19", IL20="IL-20", IL22="IL-22",
                IL24="IL-24", IL26="IL-26", IFNL3="IFN-\u03bb3")
CYT_LEVELS <- unname(CYT_LABEL)

# shared so the two arms stay comparable if placed side by side
HEATMAP_LIM  <- c(-4, 4)
HEATMAP_COLS <- c("#2166AC", "white", "#B2182B")
DIR_COLS     <- c(Activation = "#B2182B", Repression = "#2166AC")

theme_manuscript <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(plot.title   = element_text(size = base_size + 1),
          plot.caption = element_text(size = base_size - 3, colour = "grey40", hjust = 0),
          plot.caption.position = "plot",
          panel.grid.minor = element_blank())
}

# retries: SMB writes fail intermittently
safe_write <- function(x, path, n = 3) {
  for (i in seq_len(n)) {
    ok <- tryCatch({ write.csv(x, path, row.names = FALSE); TRUE },
                   error = function(e) { Sys.sleep(2); FALSE })
    if (ok) return(invisible(TRUE))
  }
  warning("could not write ", path); invisible(FALSE)
}

svg_dev <- function(f, width, height)
  svglite::svglite(f, width = width, height = height, fix_text_size = FALSE)

save_gg <- function(plot, name, width = 8, height = 6, subdir = NULL) {
  d <- if (is.null(subdir)) FIG_DIR else file.path(FIG_DIR, subdir)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  for (ext in c("png","pdf","svg")) {
    f <- file.path(d, paste0(name, ".", ext))
    dev <- switch(ext, png = "png", pdf = grDevices::pdf,
                  svg = function(filename, width, height, ...)
                          svg_dev(filename, width, height))
    try(ggsave(f, plot, width = width, height = height, dpi = 300, device = dev),
        silent = TRUE)
    sz <- file.size(f)
    if (is.na(sz) || sz < 1000) warning(sprintf("%s looks empty", basename(f)))
  }
  invisible(plot)
}

# also saves the object so panels can be re-laid-out without rerunning
save_panel <- function(plot, name, width = 8, height = 6, subdir = NULL) {
  save_gg(plot, name, width, height, subdir)
  dir.create(file.path(STATE_DIR, "panels"), showWarnings = FALSE, recursive = TRUE)
  saveRDS(plot, file.path(STATE_DIR, "panels", paste0(name, ".rds")))
  invisible(plot)
}

message(sprintf("setup: %s (%s) -> %s", DATASET, CFG$label, OUT_DIR))
GEOF

echo "updated. staged changes:"
git add -A
git status --short
