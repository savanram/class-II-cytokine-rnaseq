# Run once. Only installs what's missing.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

cran <- c("tidyverse", "ggrepel", "cowplot", "conflicted", "svglite", "config",
          "here", "rmarkdown", "patchwork", "msigdbr", "TOSTER",
          "circlize", "ggpubr")

bioc <- c("edgeR", "limma", "AnnotationDbi", "org.Hs.eg.db",
          "clusterProfiler", "enrichplot", "ComplexHeatmap",
          "progeny", "dorothea", "decoupleR")

miss <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]

if (length(miss(cran))) install.packages(miss(cran), repos = "https://cloud.r-project.org")
if (length(miss(bioc))) BiocManager::install(miss(bioc), ask = FALSE, update = FALSE)

# these four link against the same C libs; a partial upgrade breaks svglite
# silently. reinstall together from source if svg output comes out empty.
# install.packages(c("systemfonts","textshaping","ragg","svglite"), type = "source")

pkgs <- c(cran, bioc)
still <- miss(pkgs)
if (length(still)) warning("missing: ", paste(still, collapse = ", "))

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

print(data.frame(
  package = pkgs,
  version = vapply(pkgs, function(p)
    if (requireNamespace(p, quietly = TRUE))
      as.character(packageVersion(p)) else NA_character_, character(1)),
  row.names = NULL))
