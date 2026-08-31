# Full rebuild. 01-05 run per dataset; 06-08 run once across both.

library(rmarkdown)

per_dataset <- sprintf("analysis/%s.Rmd", c(
  "01_wrangle_counts", "02_deg", "03_heatmap_all_contrasts",
  "04_clustering", "05_tf_enrichment_lollipop"))

shared <- sprintf("analysis/%s.Rmd", c(
  "06_progeny_bubble", "07_cross_dataset", "08_assemble_figures"))

for (ds in c("P560-2", "P560-4")) {
  for (nb in per_dataset) {
    message("=== ", ds, " :: ", basename(nb), " ===")
    render(nb, params = list(dataset = ds),
           output_file = sprintf("%s_%s.html",
                                 tools::file_path_sans_ext(basename(nb)), ds),
           output_dir = "docs", envir = new.env())
  }
}

for (nb in shared) {
  message("=== shared :: ", basename(nb), " ===")
  render(nb, output_dir = "docs", envir = new.env())
}
