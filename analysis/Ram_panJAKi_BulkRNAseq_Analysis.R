# ---- Rscript for performing RNAseq analysis of cytokine stimulated HaCat cells +/- pan-Janus Kinase (panJAKi) inhibitor ----
# Author: Isaac M. Barber-Axthelm, DVM, MS, PhD
# Institution: University of Washington 
# Date: August 27, 2026 

# Note: # IFNb dominates the transcriptional response and is excluded from the downstream analysis to better resolve effects of IL cytokines and IFNl3.

# Study design: 
#   - Cell type: HaCaT cells
#   - Pretreatment conditions: Mock, panJAKi
#   - Stimulation conditions: Mock, IL-19, IL-20, IL-22, IL-24, IL-26, IFNl3
#   - Number of replicates per condition: 3

# This script performs the following analyses visualizations: 
#   1) Initial data processing of count matrices and sample annotations.
#   2) Aligning ensemblID's to HGNC symbols and associated gene information. 
#   3) DE analysis following the edgeR glmQLF pipeline between control and cytokine treated conditions, binned by pre-treatment groups [Mock, panJAKi]).
#       BCV and QL QC plots are generated. 
#   4) Hierarchical clustering of significant deferentially expressed genes (DEGs) across contrast groups, and plotting results as a heatmap.
#   5) Gene cluster-based transcription factor enrichment analysis using the DoRothEA viper reference gene set.
#       Significant results are plotted as lollipop plots per cluster, and top TFs are incorporated into the heatmap.
#   6) Pathway enrichment analysis comparisons using sample level pathway activity scores derived from PROGENy, between pre-treatment conditions, binned by stimulation condition [Mock, panJAKi]). 
#       Significant results are plotted as bar plots across stimulation groups, and dot plots of pathways of interest (JAK-STAT, MAPK, NF-kB). 

# Required files to run this script: 
#   1) P560-4_gene_counts.csv [supplied by authors] 
#   3) P560-4_sample_metadata.csv [supplied by authors] 
#   4) P560-4_annotations.csv [supplied by authors] 
#   5) BEST_viperRegulon.rdata [Source from DoRothEA GitHub]
#     - "https://github.com/deeenes/DoRothEA/blob/master/data/TFregulons/Robjects_VIPERformat/consensus/BEST_viperRegulon.rdata 

# ---- Load necessary packages ----
pkgs <- c("magrittr", "tidyverse", "utils", "tidyr", "dplyr", "purrr", "stringr", "tibble", "here", "stats", "grid", "scales", "rstatix", "readxl", 
          "biomaRt", "AnnotationHub", "AnnotationDbi", 
          "edgeR", "limma", "clusterProfiler", 
          "ggplot2", "ggdendro", "dendextend", "ComplexHeatmap", "circlize", "grDevices"
          ) 
invisible(lapply(pkgs, library, character.only = TRUE)) 

# ---- Define study code (from sequencing facility) ---- 
study.code <- "P560-4" 

# ---- Define Color annotation --- 
row_cluster_colors <- c("1" = "#F8766D", "2" = "#7CAE00", "3" = "#00BFC4", "4" = "#C77CFF", "5" = "#E68613", "6" = "#0CB702", "7" = "#ED68ED", "8" = "#00B8E7")
cytokine_colors <- c("IL-19" = "magenta3",  "IL-20" = "cyan3", "IL-22" = "purple", "IL-24" = "forestgreen", "IL-26" = "#8494FF", "IFNl3" = "red4") 
pretreatment_colors <- c("panJAKi" = "black", "Mock" = "#ededed") 

# ---- Define number of gene clusters for cutree----
k_row <- 3 

if (!k_row %% 1 == 0) {
  stop(paste0("Non-integer value selected for cutree. k_row value: ", k_row))
} else {
  message(paste0("Number of gene clusters for cutree: ", k_row))  
} 

# ---- Define significance thresholds ---- 
p.value.threshold <- 0.05
logFC.threshold <- log2(2) 
glm.treat.logFC.threshold <- log2(1.5) 
count.threshold <- 3 

# ---- Define plot parameters ---- 
# Plot DPI
plot.dpi <- 600 

# Bar plots
bar.height <- 9 
bar.width <- 12

# Dot plots
dot.height <- 8 
dot.width <- 14 
mad.limit <- 4 # Color scaling cutoff

# Heatmaps
heatmap.height <- 14
heatmap.width <- 12 
heatmap.legend.position <- "bottom" 
heatmap.legend.grouping <- "original" 

# Lollipop plots 
lollipop.height <- 10 
lollipop.width <- 14 

# DE analysis QC plots
qc.plot.height <- 8 
qc.plot.width <- 10 

# ---- Define the cytokine order for plotting ---- 
cytokine_order <- c("IL-19", "IL-20", "IL-22", "IL-24", "IL-26", "IFNl3") 

# ---- Load the PROGENy reference model ---- 
# Specify the number of genes to return per pathway
n.genes <- 100

# Species 
progeny.species <- "human"

message(paste0("Retreiving the ", progeny.species, " model for PROGENy pathway enrichment analysis;\n", 
               n.genes, "retreived per pathway.")) 

# Retrieve the model 
progeny_net <- decoupleR::get_progeny(organism = progeny.species, top = n.genes) 

# ---- Establish a connection with the human ensembl biomart database ----
query_biomart <- function(attributes, filter, values, dataset = "hsapiens_gene_ensembl") {
  mirrors <- c("www", "useast", "uswest", "asia")
  
  for (mirror in mirrors) {
    message("Trying mirror: ", mirror)
    tryCatch({
      mart <- biomaRt::useEnsembl(
        biomart = "ensembl",
        dataset = dataset,
        mirror  = mirror 
      ) 
      
      attribute.reference <- biomaRt::listAttributes(mart = mart) %>% dplyr::pull(name) 
      
      if(any(!attributes %in% attribute.reference)) {
        wrong_attributes <- setdiff(attributes, attribute.reference) 
        stop(paste0("The following attributes are not available in the biomart reference list: ", 
                    paste(wrong_attributes, collapse = ", ")))
      }
      
      message("Connected via: ", mirror, ". Running query...")
      result <- biomaRt::getBM(
        attributes = attributes,
        filter     = filter,
        values     = values,
        mart       = mart
      )
      message("Query succeeded via: ", mirror)
      return(list(result = result, mart = mart))
      
    }, error = function(e) {
      message("Mirror '", mirror, "' failed: ", conditionMessage(e))
    })
  }
  stop("All BioMart mirrors failed. Try again later.")
} 

# ---- Set project paths and load count/metadata -----
# Portable RNA-seq input paths for a reproducible project
# Expected project layout:
#   project_root/
#   └── P560-4/
#       ├── data/
#       │   ├── P560-4_gene_counts.csv 
#       │   ├── P560-4_annotations.csv
#       │   └── P560-4_sample_metadata.csv
#       └── results/ 

# ---- Establish project root: prefer explicit file-based detection, fall back to heuristic detection ----
this_script <- "Ram_panJAKi_BulkRNAseq_Analysis.R"

root_method <- tryCatch({
  here::i_am(this_script)
  "explicit (i_am)"
}, error = function(e) {
  warning(
    "Could not establish project root via i_am() using '", this_script, "'.\n",
    "Falling back to heuristic root detection (here::here()), which looks for the ",
    "nearest .git, .Rproj, or .here marker above the working directory.\n",
    "This fallback is less reliable if this script sits inside multiple nested ",
    "project copies, since it does not verify the script's own location.\n",
    "Original error: ", conditionMessage(e),
    call. = FALSE
  )
  "heuristic (here() fallback)"
})

project_root <- here::here()

if (!dir.exists(project_root)) {
  stop("Resolved project root does not exist on disk: ", project_root, call. = FALSE)
}

# Sanity check 
if (root_method == "heuristic (here() fallback)") {
  expected_script_path <- file.path(project_root, this_script)
  if (!file.exists(expected_script_path)) {
    warning(
      "Fallback root does not contain the expected script at: ", expected_script_path, "\n",
      "This suggests the resolved root may be the wrong project directory. ",
      "Proceed with caution, or set the working directory manually before re-running.",
      call. = FALSE
    )
  }
}

message(sprintf("Project root resolved to: %s  [method: %s]", project_root, root_method))

# Optional override for users with data stored outside the repository:
#   RNASEQ_DATA_DIR=/path/to/data Rscript code/RNAseq/pathway_analysis.R
# The default is project-relative and portable.
data_root <- Sys.getenv(
  "RNASEQ_DATA_DIR",
  unset = file.path(project_root, study.code, "data")
) 

# Create the expected input and output directories recursively.
dir.create(data_root, recursive = TRUE, showWarnings = FALSE) 

# Results directory 
working_directory <- file.path(project_root, study.code, "results")
dir.create(working_directory, recursive = TRUE, showWarnings = FALSE) 

# DoRothEA Gene set reference directory 
DoRothEA_directory <- file.path(project_root, study.code, "DoRothEA") 
dir.create(DoRothEA_directory, recursive = TRUE, showWarnings = FALSE) 

input_paths <- c(
  counts = file.path(data_root, paste0(study.code, "_gene_counts.csv")),
  metadata = file.path(data_root, paste0(study.code, "_sample_metadata.csv")),
  annotation = file.path(data_root, paste0(study.code, "_annotations.csv")), 
  DoRothEA_reference = file.path(DoRothEA_directory, "BEST_viperRegulon.rdata")  
)

check_input_paths <- function(paths) {
  paths[!file.exists(paths)]
}

missing_paths <- check_input_paths(input_paths)

if (length(missing_paths) > 0L) {
  message(
    "Required file(s) are missing.\n",
    "Add the files to the following project-relative directory/directories:\n",
    paste(
      sprintf("- %s: %s", names(missing_paths), missing_paths),
      collapse = "\n"
    ),
    "\n\nExpected data directory: ", normalizePath(
      data_root,
      winslash = "/",
      mustWork = FALSE
    )
  ) 
  
  # If the DoRothEA viper-regulon reference is among the missing files,
  # point the user to the public repository it can be downloaded from.
  if ("DoRothEA_reference" %in% names(missing_paths)) {
    message(
      "\nThe DoRothEA viper-regulon reference object was not found.\n",
      "Download 'BEST_viperRegulon.rdata' from the public repository:\n",
      "https://github.com/deeenes/DoRothEA/blob/master/data/TFregulons/Robjects_VIPERformat/consensus/BEST_viperRegulon.rdata\n",
      "and place it at: ", normalizePath(
        input_paths[["DoRothEA_reference"]],
        winslash = "/",
        mustWork = FALSE
      )
    )
  } 
  
  # In an interactive R/RStudio session, allow the user to add the files
  # and then recheck before the analysis proceeds.
  if (interactive()) {
    invisible(readline(
      "Add the missing CSV files, then press <Enter> to recheck: "
    ))
    missing_paths <- check_input_paths(input_paths)
  }
  
  # Stop before analysis if files are still unavailable. This is also the
  # behavior for noninteractive Rscript execution.
  if (length(missing_paths) > 0L) {
    stop(
      "Analysis stopped because required input file(s) are missing:\n",
      paste(
        sprintf("- %s: %s", names(missing_paths), missing_paths),
        collapse = "\n"
      ),
      "\n\nPlace the files in the documented data/ directory and rerun the script."
    )
  }
}

# All required inputs are present; load them.
counts_raw <- utils::read.csv(
  file = input_paths[["counts"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_metadata <- utils::read.csv(
  file = input_paths[["metadata"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)

annotation <- utils::read.csv(
  file = input_paths[["annotation"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
) 

load(file = input_paths[["DoRothEA_reference"]])  
if (!exists("viper_regulon")) {
  stop("Loaded RData did not contain an object named 'viper_regulon' — check the reference file.")
} 

# ---- Join sample metadata and library annotations ---- 
sample_metadata <- sample_metadata %>%
  dplyr::left_join(
    annotation %>%
      dplyr::mutate(sample_name = sub("^[0-9]+_", "", sample_name)),
    by = c("Sample.name" = "sample_name")
  ) 

# ---- Generate a gene metadata object ----
gene_metadata <- counts_raw %>% dplyr::select(geneName, gene_symbol)

# Prepare counts object for DE analysis 
counts_raw <- counts_raw %>%
  tibble::column_to_rownames(var = "geneName") %>% 
  dplyr::select(-gene_symbol) 
colnames(counts_raw) <- sub("_.*", "", colnames(counts_raw))

# ---- Query biomart and AnnotationHub to generate row metadata ----
attrs <- c("ensembl_gene_id", "description", "gene_biotype", "entrezgene_id", "hgnc_symbol")

biomart_out <- query_biomart(
  attributes = attrs,
  filter     = "ensembl_gene_id",
  values     = rownames(counts_raw)
)

row_metadata_biomart <- biomart_out$result 

gene_metadata <- gene_metadata %>%
  dplyr::left_join(row_metadata_biomart, by = c("geneName" = "ensembl_gene_id")) %>% 
  dplyr::mutate(
    gene_symbol = dplyr::coalesce(gene_symbol, hgnc_symbol)
  ) %>% 
  dplyr::select(-hgnc_symbol)

missing_gene_annotation <- gene_metadata %>%
  dplyr::filter(is.na(gene_symbol)) %>%
  dplyr::pull(geneName)

ah <- AnnotationHub::AnnotationHub()
edb_query <- AnnotationHub::query(ah, c("EnsDb", "Homo sapiens"))

ah_ids    <- names(edb_query)
ah_titles <- edb_query$title
message("Sweeping across ", length(ah_ids), " Ensembl versions...")

sweep_results <- list()

for (i in seq_along(ah_ids)) {
  ah_id <- ah_ids[i]
  message("Trying ", ah_titles[i], " (", ah_id, ")...")
  
  tryCatch({
    edb_i <- edb_query[[ah_id]]
    
    result_i <- AnnotationDbi::select(
      edb_i,
      keys    = missing_gene_annotation,
      keytype = "GENEID",
      columns = c("GENEID", "SYMBOL", "GENEBIOTYPE", "ENTREZID", "DESCRIPTION")
    ) %>%
      dplyr::filter(!is.na(GENEBIOTYPE)) %>%
      dplyr::mutate(source_version = ah_titles[i])
    
    if (nrow(result_i) > 0) {
      message("  -> Recovered ", nrow(result_i), " annotations")
      sweep_results[[ah_id]] <- result_i
    }
    
  }, error = function(e) {
    message("  -> Failed: ", conditionMessage(e))
  })
}

ah_sweep_combined <- dplyr::bind_rows(sweep_results) %>%
  dplyr::group_by(GENEID) %>%
  dplyr::slice_tail(n = 1) %>%
  dplyr::ungroup()

message("Total unique genes recovered: ", dplyr::n_distinct(ah_sweep_combined$GENEID),
        " of ", length(missing_gene_annotation))

ah_result_clean <- ah_sweep_combined %>%
  dplyr::rename(
    geneName      = GENEID,
    gene_symbol   = SYMBOL,
    gene_biotype  = GENEBIOTYPE,
    entrezgene_id = ENTREZID,
    description   = DESCRIPTION
  ) %>%
  dplyr::select(-source_version)

gene_metadata <- gene_metadata %>%
  dplyr::rows_patch(ah_result_clean, by = "geneName", unmatched = "ignore")

still_missing <- gene_metadata %>% dplyr::filter(is.na(gene_symbol))

if (nrow(still_missing) > 0) {
  message("Still unresolved after sweep: ", nrow(still_missing))
} else {
  message("All ensembl IDs successfully mapped.")
}

gene_metadata <- gene_metadata %>%
  dplyr::distinct(geneName, .keep_all = TRUE)

# ---- Verify column alignment ----
stopifnot(all(colnames(counts_raw) == sample_metadata$libid))

# ---- Parse metadata into structured factors ----
sample_metadata <- sample_metadata %>%
  dplyr::mutate(
    cytokine = dplyr::case_when(
      grepl("IL-19",       Condition) ~ "IL-19",
      grepl("IL-20",       Condition) ~ "IL-20",
      grepl("IL-22",       Condition) ~ "IL-22",
      grepl("IL-24",       Condition) ~ "IL-24",
      grepl("IL-26",       Condition) ~ "IL-26",
      grepl("IFN\u03bb3", Condition) ~ "IFNl3", 
      grepl("IFN\u03b2", Condition) ~ "IFNb", 
      TRUE ~ "Mock"
    ),
    arm = dplyr::case_when(
      grepl("Jak inhibitor", Condition) ~ "JakI",
      grepl("stim",          Condition) ~ "Stim",
      TRUE ~ "Mock"
    ),
    Replicate = factor(Replicate)
  ) %>% 
  # Build a single condition factor encoding every unique group 
  dplyr::mutate(
    condition_level = dplyr::case_when(
      arm == "Mock" ~ "Mock", 
      TRUE          ~ paste0(cytokine, "_", arm)
    ), 
    condition_level = stringr::str_replace_all(condition_level, "-", ""), 
    condition_level = factor(condition_level)
  ) 

# ---- Exclude IFNb samples ----
message("Excluding IFNb samples from analysis...")
message("  Samples before exclusion: ", nrow(sample_metadata))

ifnb_libids <- sample_metadata %>%
  dplyr::filter(cytokine == "IFNb") %>%
  dplyr::pull(libid)

message("  IFNb samples to exclude: ", paste(ifnb_libids, collapse = ", "))

sample_metadata <- sample_metadata %>%
  dplyr::filter(cytokine != "IFNb") %>%
  dplyr::mutate(condition_level = droplevels(condition_level))

counts_raw <- counts_raw[, sample_metadata$libid, drop = FALSE]

message("  Samples after exclusion: ", nrow(sample_metadata))
message("  Remaining cytokines: ",
        paste(sort(unique(sample_metadata$cytokine)), collapse = ", "))

stopifnot(all(colnames(counts_raw) == sample_metadata$libid))
message("  PASS: counts_raw columns realigned after IFNb exclusion.")

print(sample_metadata %>% dplyr::select(libid, Sample.name, cytokine, arm, Replicate, condition_level))

message("Condition levels:")
print(levels(sample_metadata$condition_level)) 

# ---- Sample Assignment and Contrast Validation Checks ----
# ---- Check: Verify every sample has a valid cytokine, arm, and condition_level ----
message("=== Check: Metadata completeness ===")

metadata_check <- sample_metadata %>%
  dplyr::summarise(
    n_samples          = dplyr::n(),
    n_missing_cytokine = sum(is.na(cytokine)),
    n_missing_arm      = sum(is.na(arm)),
    n_missing_condition = sum(is.na(condition_level)),
    n_unexpected_cytokine = sum(!cytokine %in% c("IL-19", "IL-20", "IL-22", "IL-24",
                                                 "IL-26", "IFNl3", "Mock")), 
    n_unexpected_arm      = sum(!arm %in% c("Stim", "JakI", "Mock"))
  )
print(metadata_check)

if(any(metadata_check[, -1] > 0)) {
  stop("Metadata completeness check failed — review sample_metadata parsing above.")
} else {
  message("  PASS: All samples have valid cytokine, arm, and condition_level assignments.")
}

# ---- Check: Verify expected sample counts per condition ----
message("=== Check: Sample counts per condition ===")

condition_counts <- sample_metadata %>%
  dplyr::count(cytokine, arm, Replicate, condition_level) %>%
  dplyr::arrange(cytokine, arm, Replicate)
print(condition_counts)

# Each condition should have exactly 1 sample (3 replicates per cytokine×arm)
if(any(condition_counts$n != 1)) {
  stop("Duplicate samples detected in the same condition×replicate cell — check metadata.")
} else {
  message("  PASS: Each condition × replicate has exactly 1 sample.")
}

# Each non-Mock condition should have exactly 3 replicates
replicate_counts <- sample_metadata %>%
  dplyr::count(cytokine, arm) %>%
  dplyr::arrange(cytokine, arm)
print(replicate_counts)

expected_n <- 3
non_mock_check <- replicate_counts %>%
  dplyr::filter(!(cytokine == "Mock" & arm == "Mock")) %>%
  dplyr::filter(n != expected_n)

mock_check <- replicate_counts %>%
  dplyr::filter(cytokine == "Mock" & arm == "Mock") %>%
  dplyr::filter(n != expected_n)

if(nrow(non_mock_check) > 0) {
  message("  WARNING: Some non-Mock conditions do not have exactly 3 replicates:")
  print(non_mock_check)
} else {
  message("  PASS: All non-Mock conditions have exactly 3 replicates.")
}

if(nrow(mock_check) > 0) {
  message("  WARNING: Mock condition does not have exactly 3 replicates.")
} else {
  message("  PASS: Mock has exactly 3 replicates.")
}

# ---- Check: Verify libid-to-condition alignment between metadata and counts ----
message("=== Check: Column order alignment between counts and metadata ===")

# These must be identical in both order and content
counts_cols   <- colnames(counts_raw)
metadata_libs <- sample_metadata$libid

if(!identical(counts_cols, metadata_libs)) {
  mismatches <- data.frame(
    position   = seq_along(counts_cols),
    counts_col = counts_cols,
    metadata_libid = metadata_libs,
    match = counts_cols == metadata_libs
  )
  print(mismatches[!mismatches$match, ])
  stop("Column order mismatch between counts_raw and sample_metadata — DE analysis would use wrong groups.")
} else {
  message("  PASS: counts_raw columns exactly match sample_metadata libid order.")
}

# ---- Check: Verify condition_level matches the libid naming convention ----
message("=== Check: Sample name to condition cross-reference ===")

# Print a cross-reference table for manual inspection
cross_ref <- sample_metadata %>%
  dplyr::select(libid, Sample.name, cytokine, arm, Replicate, condition_level) %>%
  dplyr::arrange(cytokine, arm, Replicate)
print(as.data.frame(cross_ref))

# Flag any Mock samples that were assigned a cytokine
mock_cytokine_mismatch <- sample_metadata %>%
  dplyr::filter(arm == "Mock" & cytokine != "Mock")
if(nrow(mock_cytokine_mismatch) > 0) {
  message("  WARNING: Mock arm samples assigned non-Mock cytokine:")
  print(mock_cytokine_mismatch)
} else {
  message("  PASS: All Mock arm samples have cytokine = Mock.")
}

# Flag any stimulated/JakI samples assigned cytokine = Mock
stim_mock_mismatch <- sample_metadata %>%
  dplyr::filter(arm != "Mock" & cytokine == "Mock")
if(nrow(stim_mock_mismatch) > 0) {
  message("  WARNING: Non-Mock arm samples assigned cytokine = Mock:")
  print(stim_mock_mismatch)
} else {
  message("  PASS: All non-Mock arm samples have a valid cytokine assignment.")
} 

# ---- Build global design matrix (no intercept) ----
# ~ 0 + condition gives one coefficient per condition,
# making explicit contrasts straightforward to specify.
design_global <- model.matrix(~ 0 + condition_level + Replicate, data = sample_metadata)
colnames(design_global) <- c(levels(sample_metadata$condition_level),
                             paste0("Replicate", levels(sample_metadata$Replicate)[-1])) 

message("Design matrix columns: ", paste(colnames(design_global), collapse = ", "))
message("Design rank: ", qr(design_global)$rank)
message("Residual df: ", nrow(design_global) - qr(design_global)$rank) 

# ---- Check: Verify design matrix construction ----
message("=== Check: Design matrix validation ===")

# Every condition_level should map to exactly one column in the design
expected_replicate_cols <- paste0("Replicate", levels(sample_metadata$Replicate)[-1])
expected_all <- c(levels(sample_metadata$condition_level), expected_replicate_cols)

if(!identical(sort(colnames(design_global)), sort(expected_all))) {
  stop("Design matrix mismatch.")
} else {
  message("  PASS: Design matrix columns match all expected condition levels.")
}

# Check each sample maps to exactly one design column
condition_cols <- colnames(design_global) %in% levels(sample_metadata$condition_level)
design_row_sums_condition <- rowSums(design_global[, condition_cols, drop = FALSE])
if(any(design_row_sums_condition != 1)) {
  stop("Design matrix row sum check failed (condition-level columns).")
} else {
  message("  PASS: Every sample maps to exactly one condition column.")
}

# Verify each condition has the expected number of samples
design_col_sums <- colSums(design_global[, condition_cols, drop = FALSE]) 
message("  Samples per condition in design matrix:")
print(sort(design_col_sums))

expected_per_condition <- 3
wrong_counts <- design_col_sums[design_col_sums != expected_per_condition]
if(length(wrong_counts) > 0) {
  message("  WARNING: These conditions do not have exactly ", expected_per_condition, " samples:")
  print(wrong_counts)
} else {
  message("  PASS: Every condition has exactly ", expected_per_condition, " samples.")
}

# ---- Check: Verify contrast construction ----
message("=== Check: Contrast vector validation ===")

# Build and inspect every contrast before running DE
contrast_audit <- tidyr::expand_grid(
  cytokine_id = setdiff(unique(sample_metadata$cytokine), "Mock"),
  test_arm    = c("Stim", "JakI")
) %>%
  dplyr::mutate(
    clean_cytokine = stringr::str_replace_all(cytokine_id, "-", ""),
    condition_name = paste0(clean_cytokine, "_", test_arm),
    contrast_label = paste0(clean_cytokine, "_", test_arm, "_vs_Mock"),
    condition_in_design = condition_name %in% colnames(design_global),
    n_samples_in_condition = colSums(design_global)[condition_name],
    mock_in_design = "Mock" %in% colnames(design_global) 
  )

print(contrast_audit, n = Inf)

# Flag any intended contrasts where the condition is missing from the design
missing_conditions <- contrast_audit %>% dplyr::filter(!condition_in_design)
if(nrow(missing_conditions) > 0) {
  message("  WARNING: The following contrasts reference conditions not in the design:")
  print(missing_conditions %>% dplyr::select(cytokine_id, test_arm, condition_name))
} else {
  message("  PASS: All intended contrasts reference conditions present in the design.")
} 

# ---- Restrict analysis to protein-coding genes ----
protein_coding_genes <- gene_metadata %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  dplyr::pull(geneName)

message(sprintf("Protein-coding genes in gene_metadata: %d of %d total",
                length(protein_coding_genes), nrow(gene_metadata)))

# Subset counts_raw to protein-coding genes present in counts_raw
protein_coding_in_counts <- intersect(protein_coding_genes, rownames(counts_raw))

message(sprintf("Protein-coding genes present in counts_raw: %d",
                length(protein_coding_in_counts)))

counts_raw_pc <- counts_raw[protein_coding_in_counts, , drop = FALSE]

message(sprintf("counts_raw reduced from %d to %d genes (protein-coding only)",
                nrow(counts_raw), nrow(counts_raw_pc))) 

# ---- DGEList, filter, normalize across ALL samples ----
y_global <- edgeR::DGEList(counts = counts_raw_pc, group = sample_metadata$condition_level)
keep <- edgeR::filterByExpr(y_global, design = design_global)
message(sprintf("filterByExpr: kept %d / %d genes (%.1f%%)",
                sum(keep), nrow(y_global), 100 * mean(keep))) 
y_global  <- y_global[keep, , keep.lib.sizes = FALSE]
y_global  <- edgeR::calcNormFactors(y_global, method = "TMM") 

# ---- Check: Library size inspection ----
message("=== Check: Library size inspection ===")

lib_sizes <- data.frame(
  libid         = colnames(counts_raw_pc),
  condition     = sample_metadata$condition_level,
  library_size  = colSums(counts_raw_pc)
) %>%
  dplyr::arrange(library_size)

message("  Library size summary:")
print(summary(lib_sizes$library_size))

# Flag samples with library sizes more than 2 SD from the mean
lib_mean <- mean(lib_sizes$library_size)
lib_sd   <- sd(lib_sizes$library_size)
outlier_libs <- lib_sizes %>%
  dplyr::filter(abs(library_size - lib_mean) > 2 * lib_sd)

if(nrow(outlier_libs) > 0) {
  message("  WARNING: The following samples have outlier library sizes (>2 SD from mean):")
  print(outlier_libs)
} else {
  message("  PASS: No outlier library sizes detected.")
}

# ---- Check: Normalization factor inspection ----
message("=== Check: TMM normalization factor inspection ===")

norm_factors <- data.frame(
  libid        = colnames(y_global$counts),
  condition    = sample_metadata$condition_level,
  norm_factor  = y_global$samples$norm.factors
) %>%
  dplyr::arrange(norm_factor)

message("  Normalization factor summary:")
print(summary(norm_factors$norm_factor))

# Flag extreme normalization factors (outside 0.5 – 2.0 is unusual)
extreme_nf <- norm_factors %>%
  dplyr::filter(norm_factor < 0.5 | norm_factor > 2.0)

if(nrow(extreme_nf) > 0) {
  message("  WARNING: Extreme normalization factors detected — these samples may be outliers:")
  print(extreme_nf)
} else {
  message("  PASS: All normalization factors within normal range (0.5 – 2.0).")
}

# ---- Global dispersion estimation ----
# Estimating dispersion across all samples. 
# available to the empirical Bayes shrinkage procedure (Robinson & Smyth 2007).
y_global <- edgeR::estimateDisp(y_global, design_global, robust = TRUE)
message(sprintf("Global common dispersion: %.4f (BCV: %.4f)",
                y_global$common.dispersion,
                sqrt(y_global$common.dispersion)))

# BCV plot — inspect dispersion-mean trend across all samples
# Using a ggplot2-based plotting method

# Extract data from the DGEList object
bcv_df <- data.frame(
  AveLogCPM = y_global$AveLogCPM,
  tagwise_BCV = sqrt(y_global$tagwise.dispersion),
  trended_BCV = sqrt(y_global$trended.dispersion),
  common_BCV = sqrt(y_global$common.dispersion)
)

# Order by AveLogCPM for a smooth trend line
bcv_df <- bcv_df[order(bcv_df$AveLogCPM), ]

# Common BCV value (single scalar for horizontal line)
common_bcv_value <- sqrt(y_global$common.dispersion)

bcv_plot <- ggplot2::ggplot(bcv_df, ggplot2::aes(x = AveLogCPM)) +
  # Tagwise dispersions (individual gene points)
  ggplot2::geom_point(
    ggplot2::aes(y = tagwise_BCV),
    color  = "black",
    alpha  = 0.2,
    size   = 0.4,
    shape  = 16
  ) +
  # Trended dispersion line
  ggplot2::geom_line(
    ggplot2::aes(y = trended_BCV),
    color     = "#2166AC",
    linewidth = 0.8
  ) +
  # Common dispersion horizontal line
  ggplot2::geom_hline(
    yintercept = common_bcv_value,
    color      = "#D6604D",
    linewidth  = 0.8,
    linetype   = "solid"
  ) +
  ggplot2::scale_x_continuous(
    name   = expression("Average log"[2]*" CPM"),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    name   = "Biological coefficient of variation",
    expand = ggplot2::expansion(mult = c(0.02, 0.05))
  ) +
  ggplot2::annotate(
    "text",
    x     = Inf, # Right border
    y     = Inf, # Top border
    label = sprintf("Common BCV = %.3f", common_bcv_value),
    hjust = 1.1, # Small inset from right edge
    vjust = 1.5, # Small inset from top edge
    size  = 3,
    color = "#D6604D"
  ) +
  ggplot2::labs(title = sprintf("BCV plot: global dispersion estimate (%d samples)",
                                nrow(sample_metadata)))  +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title   = ggplot2::element_text(size = 12, face = "bold", color = "black"),
    axis.title   = ggplot2::element_text(size = 12, face = "bold", color = "black"),
    axis.text    = ggplot2::element_text(size = 10, face = "plain", color = "black"),
    panel.grid.minor = ggplot2::element_blank()
  )

print(bcv_plot) 

# ---- Fit global QL model ----
fit_global <- edgeR::glmQLFit(y_global, design_global, robust = TRUE)

# QL dispersion plot
# Using a ggplot2-based plotting method 

ql_df <- data.frame(
  AveLogCPM = fit_global$AveLogCPM,
  var_prior = fit_global$s2.prior, 
  var_post = fit_global$s2.post 
)

ql_df <- ql_df[order(ql_df$AveLogCPM), ]

ql_plot <- ggplot2::ggplot(ql_df, ggplot2::aes(x = AveLogCPM)) + 
  # Raw (posterior) QL dispersion points 
  ggplot2::geom_point(
    ggplot2::aes(y = var_post),
    color = "black",
    alpha = 0.2,
    size  = 0.4,
    shape = 16
  ) + 
  # Prior trend line (empirical Bayes prior) 
  ggplot2::geom_line(
    ggplot2::aes(y = var_prior),
    color     = "#2166AC",
    linewidth = 0.8
  ) + 
  # Reference line at 1 (no extra-Poisson variation) 
  ggplot2::geom_hline(
    yintercept = 1,
    color      = "#D6604D",
    linewidth  = 0.8,
    linetype   = "dashed"
  ) + 
  ggplot2::scale_x_continuous(
    name   = expression("Average log"[2]*" CPM"),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) + 
  ggplot2::scale_y_continuous(
    name   = "Quarter-root mean QL dispersion", 
    trans = scales::trans_new(
      name      = "quarter_root",
      transform = function(x) x^0.25,
      inverse   = function(x) x^4
    ), 
    expand = ggplot2::expansion(mult = c(0.02, 0.05))
  ) + 
  ggplot2::labs(title = sprintf("QL dispersion: global model (%d samples)",
                                nrow(sample_metadata))) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(size = 12, face = "bold", color = "black"),
    axis.title       = ggplot2::element_text(size = 12, face = "bold", color = "black"),
    axis.text        = ggplot2::element_text(size = 10, face = "plain", color = "black"),
    panel.grid.minor = ggplot2::element_blank()
  )

print(ql_plot)

# ---- Save QC plots ----
outdir_qc <- file.path(working_directory, "edgeR_RNAseq_QC")
dir.create(outdir_qc, showWarnings = FALSE, recursive = TRUE) 

# PDF
ggplot2::ggsave(
  filename = paste("ggplot", "panJAKi", "global", "BCV", "plot.pdf", sep = "_"),
  path = outdir_qc,
  plot = bcv_plot, 
  device = cairo_pdf, 
  dpi = plot.dpi,
  height = qc.plot.height,
  width = qc.plot.width,
  units = "cm",
  bg = "transparent"
) 

ggplot2::ggsave(
  filename = paste("ggplot", "panJAKi", "global", "QL", "dispersion", "plot.pdf", sep = "_"),
  path = outdir_qc,
  plot = ql_plot, 
  device = cairo_pdf, 
  dpi = plot.dpi,
  height = qc.plot.height,
  width = qc.plot.width,
  units = "cm",
  bg = "transparent"
) 

# PNG
ggplot2::ggsave(
  filename = paste("ggplot", "panJAKi", "global", "BCV", "plot.png", sep = "_"),
  path = outdir_qc,
  plot = bcv_plot, 
  device = "png", 
  dpi = plot.dpi,
  height = qc.plot.height,
  width = qc.plot.width,
  units = "cm",
  bg = "transparent"
) 

ggplot2::ggsave(
  filename = paste("ggplot", "panJAKi", "global", "QL", "dispersion", "plot.png", sep = "_"),
  path = outdir_qc,
  plot = ql_plot, 
  device = "png", 
  dpi = plot.dpi,
  height = qc.plot.height,
  width = qc.plot.width,
  units = "cm",
  bg = "transparent"
) 

# SVG
ggplot2::ggsave(
  filename = paste("ggplot", "panJAKi", "global", "BCV", "plot.svg", sep = "_"),
  path = outdir_qc,
  plot = bcv_plot, 
  device = "svg", 
  dpi = plot.dpi,
  height = qc.plot.height,
  width = qc.plot.width,
  units = "cm",
  bg = "transparent"
) 

ggplot2::ggsave(
  filename = paste("ggplot", "panJAKi", "global", "QL", "dispersion", "plot.svg", sep = "_"),
  path = outdir_qc,
  plot = ql_plot, 
  device = "svg", 
  dpi = plot.dpi,
  height = qc.plot.height,
  width = qc.plot.width,
  units = "cm",
  bg = "transparent"
) 

# ---- Define and extract contrasts of interest ----
# Each contrast compares one cytokine condition vs the three Mock replicates. 
cytokines_to_test <- setdiff(unique(sample_metadata$cytokine), "Mock")
arms_to_test      <- c("Stim", "JakI") 
stopifnot(all(cytokines_to_test %in% sample_metadata$cytokine), all(arms_to_test %in% sample_metadata$arm)) 

extract_contrast <- function(fit, design, cytokine_id, test_arm,
                             treat_lfc = glm.treat.logFC.threshold,
                             sig_fdr   = p.value.threshold,
                             verbose   = TRUE) {
  
  .say <- function(...) if (verbose) message(sprintf(...)) 
  # Clean contrast names
  clean_cytokine <- stringr::str_replace_all(cytokine_id, "-", "")
  condition_name <- paste0(clean_cytokine, "_", test_arm) 
  mock_name      <- "Mock"  
  contrast_label <- paste0(clean_cytokine, "_", test_arm, "_vs_Mock")
  
  .say("=== %s ===", contrast_label)
  
  # Verify the condition exists in the design
  if (!condition_name %in% colnames(design)) {
    .say("Condition '%s' not found in design — skipping.", condition_name)
    return(list(status = "condition_not_found",
                cytokine = cytokine_id,
                test_arm = test_arm,
                group.name = contrast_label))
  } 
  
  if (!mock_name %in% colnames(design)) {
    .say("Mock condition '%s' not found in design — skipping.", mock_name)
    return(list(status = "mock_not_found", cytokine = cytokine_id,
                test_arm = test_arm, group.name = contrast_label))
  } 
  
  # Build contrast: condition vs mock (3 replicates) 
  contrast_str <- paste0(condition_name, " - ", mock_name)
  
  contrast_vec <- limma::makeContrasts(
    contrasts = contrast_str,
    levels    = design
  )
  colnames(contrast_vec) <- contrast_label
  
  # ---- glmQLFTest: no LFC threshold ----
  qlf   <- edgeR::glmQLFTest(fit, contrast = contrast_vec)
  tt_ql <- edgeR::topTags(qlf, n = Inf)$table %>%
    dplyr::mutate(
      ensemblID = rownames(.),
      cytokine  = cytokine_id,
      test_arm  = test_arm,
      contrast  = contrast_label
    ) %>%
    dplyr::relocate(ensemblID, cytokine, test_arm, contrast)
  
  # ---- glmTreat: with LFC threshold ----
  tr <- edgeR::glmTreat(fit, contrast = contrast_vec, lfc = treat_lfc)
  tt_tr <- edgeR::topTags(tr, n = Inf)$table %>%
    dplyr::mutate(
      ensemblID = rownames(.),
      cytokine  = cytokine_id,
      test_arm  = test_arm,
      contrast  = contrast_label
    ) %>%
    dplyr::relocate(ensemblID, cytokine, test_arm, contrast)
  
  .say("  QL sig: %d | Treat sig: %d (FDR < %.2g, lfc >= %.2fx)",
       sum(tt_ql$FDR < sig_fdr, na.rm = TRUE),
       sum(tt_tr$FDR < sig_fdr, na.rm = TRUE),
       sig_fdr, 2^treat_lfc)
  
  list(
    status      = "ok",
    cytokine    = cytokine_id,
    test_arm    = test_arm,
    group.name  = contrast_label,
    table_ql    = tt_ql,
    table_treat = tt_tr
  )
}

# ---- Run all contrasts ----
contrasts_to_run <- tidyr::expand_grid(
  cytokine_id = cytokines_to_test,
  test_arm    = arms_to_test
)

results <- purrr::pmap(
  contrasts_to_run,
  function(cytokine_id, test_arm) {
    extract_contrast(
      fit         = fit_global,
      design      = design_global,
      cytokine_id = cytokine_id,
      test_arm    = test_arm,
      treat_lfc   = glm.treat.logFC.threshold,
      sig_fdr     = p.value.threshold
    )
  }
)
names(results) <- purrr::map_chr(results, ~ .x$group.name)

# ---- Merge results tables ----
merged_ql <- purrr::map_dfr(results, ~ if (.x$status == "ok") .x$table_ql)
merged_tr <- purrr::map_dfr(results, ~ if (.x$status == "ok") .x$table_treat)

# ---- Sanity check ----
merged_ql %>%
  dplyr::group_by(cytokine, test_arm) %>%
  dplyr::summarise(
    n_tested = dplyr::n(),
    n_sig_ql = sum(FDR < p.value.threshold, na.rm = TRUE),
    .groups  = "drop"
  ) %>% 
  print()

merged_tr %>%
  dplyr::group_by(cytokine, test_arm) %>%
  dplyr::summarise(
    n_tested    = dplyr::n(),
    n_sig_treat = sum(FDR < p.value.threshold, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  print()

# ---- Map gene metadata back to results ----
merged_ql <- merged_ql %>% 
  dplyr::mutate( 
    Significance = dplyr::case_when(
      abs(logFC) >= logFC.threshold & FDR < p.value.threshold ~ "Significant", 
      TRUE ~ "Not Significant"
    ) 
  ) %>% 
  dplyr::left_join(gene_metadata, by = c("ensemblID" = "geneName")) 

merged_tr <- merged_tr %>% 
  dplyr::mutate(
    Significance = dplyr::case_when(
      FDR < p.value.threshold ~ "Significant", 
      TRUE ~ "Not Significant"
    )
  ) %>% 
  dplyr::left_join(gene_metadata, by = c("ensemblID" = "geneName")) 

# ---- Save results ----
outdir_deg <- file.path(working_directory, "edgeR_glmQLF_differential_analysis")
dir.create(outdir_deg, showWarnings = FALSE, recursive = TRUE)

utils::write.csv(x = merged_ql, file = file.path(outdir_deg, "DEG_combined_glmQLF_results.csv"),   row.names = FALSE)
utils::write.csv(x = merged_tr, file = file.path(outdir_deg, "DEG_combined_glmTreat_results.csv"), row.names = FALSE) 
utils::write.csv(x = cross_ref, file = file.path(outdir_deg, "sample_cross_reference_table.csv"), row.names = FALSE) 

# --- Make differential Heatmap dataframe ---- 
rnaseq_heatmap_df <- merged_ql %>% 
  dplyr::group_by(ensemblID) %>% 
  dplyr::filter(any(abs(logFC) >= logFC.threshold & FDR < p.value.threshold)) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(ensemblID, logFC, contrast) %>% 
  tidyr::pivot_wider(names_from = contrast, values_from = logFC, values_fill = 0) %>% 
  tibble::column_to_rownames(var = "ensemblID") 

# ---- Hierarchial cluster rows (genes) ---- 
# z-score the log2FC data
df <- t(scale(t(rnaseq_heatmap_df))) 

# Compute Euclidean distances
dist_matrix <- stats::dist(df, method = "euclidean") 

# Perform hierarchical clustering using ward.D2
hc_row <- stats::hclust(dist_matrix, method = "ward.D2") 

# Convert to ggdendro
dendro_data <- ggdendro::dendro_data(hc_row)

# Plot the dendrogram
plot1 <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = dendro_data$segments,
                        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Hierarchical Clustering (Ward.D2)", 
                x = NULL, y = "Height") +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                 axis.ticks.x = ggplot2::element_blank()) 
print(plot1) 

# Cut tree and evaluate cut height 
clusters_row <- dendextend::cutree(tree = hc_row, k = k_row)
my_gene_col <- as.data.frame(clusters_row) 
colnames(my_gene_col) <- c("Cluster")
my_gene_col$Cluster <- as.factor(my_gene_col$Cluster)
cut_height <- sort(hc_row$height, decreasing = TRUE)[k_row - 1]

# Plot dendrogram with cut line
plot2 <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = dendro_data$segments,
                        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) + 
  ggplot2::geom_hline(yintercept = cut_height, color = "firebrick", linewidth = 0.75, linetype = "dashed") + 
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Hierarchical Clustering (Ward.D2)", 
                x = NULL, y = "Height") +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                 axis.ticks.x = ggplot2::element_blank()) 
print(plot2) 

# Annotate my_peak_cols w/ peak metadata
my_gene_col$geneName <- rownames(my_gene_col)
my_gene_col <- my_gene_col %>% 
  dplyr::left_join(gene_metadata, by = "geneName") 
rownames(my_gene_col) <- my_gene_col$geneName  

# Generate a vector of clusters to loop through, for later pathway enrichment analysis
cluster.list <- sort(unique(as.integer(as.character(my_gene_col$Cluster)))) 

# ---- Save row clustering results ---- 
outdir_cluster <- file.path(working_directory, "heatmap_gene_clustering")
dir.create(outdir_cluster, showWarnings = FALSE, recursive = TRUE)

# Save clustering assignments 
utils::write.csv(x = my_gene_col, file = file.path(outdir_cluster, "differential_RNAseq_Gene_Cluster_Assignments.csv"), row.names = FALSE)  

# ---- Hierarchical cluster columns (contrasts) ---- 
# z-score the log2FC data
df <- scale(t(rnaseq_heatmap_df)) 

# Compute Euclidean distances
dist_matrix <- stats::dist(df, method = "euclidean") 

# Perform hierarchical clustering using ward.D2
hc_col <- stats::hclust(dist_matrix, method = "ward.D2") 

# Convert to ggdendro
dendro_data <- ggdendro::dendro_data(hc_col)

# Plot the dendrogram
plot1 <- ggplot2::ggplot() +
  ggplot2::geom_segment(data = dendro_data$segments,
                        ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "Hierarchical Clustering (Ward.D2)", 
                x = NULL, y = "Height") +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), 
                 axis.ticks.x = ggplot2::element_blank()) 
print(plot1) 

# ---- Draw Heatmap with z-score normalization ---- 

#Prep data for heatmap
mat <- as.matrix(t(scale(t(rnaseq_heatmap_df)))) 

#Column Annotation 
my_sample_col <- data.frame(
  contrast = colnames(mat) 
) %>% 
  dplyr::mutate(
    Cytokine = str_extract(contrast, "^[^_]+"), 
    Pretreatment = if_else(str_detect(contrast, "^[^_]+_JakI_"), "panJAKi", "Mock") 
  ) %>% 
  dplyr::mutate(
    Cytokine = dplyr::case_when( 
      Cytokine == "IL19" ~ "IL-19", 
      Cytokine == "IL20" ~ "IL-20", 
      Cytokine == "IL22" ~ "IL-22", 
      Cytokine == "IL24" ~ "IL-24", 
      Cytokine == "IL26" ~ "IL-26", 
      TRUE               ~ Cytokine
    )
  ) %>% 
  dplyr::mutate(
    Cytokine = factor(Cytokine, levels = cytokine_order)
  )
rownames(my_sample_col) <- my_sample_col$contrast

my_sample_col <- my_sample_col[colnames(mat), , drop = FALSE] 

column_ha<- ComplexHeatmap::HeatmapAnnotation( 
  cytokine = my_sample_col$Cytokine,  
  pretreatment = my_sample_col$Pretreatment, 
  annotation_legend_param = list( 
    cytokine = list(title = "Cytokine", direction = "horizontal", 
                    title_gp = grid::gpar(fontsize = 14, 
                                          color = "black", 
                                          fontface = "bold"), 
                    title_position = "leftcenter", nrow = 1, 
                    legend_gp = grid::gpar(fontsize = 14, 
                                           color = "black", 
                                           fontface = "plain")
    ), 
    pretreatment = list(title = "Pre-Tx", direction = "horizontal", 
                        title_gp = grid::gpar(fontsize = 14, 
                                              color = "black", 
                                              fontface = "bold"), 
                        title_position = "leftcenter", nrow = 1, 
                        legend_gp = grid::gpar(fontsize = 14, 
                                               color = "black", 
                                               fontface = "plain")
    )
  ), 
  col = list( 
    cytokine = cytokine_colors, 
    pretreatment = pretreatment_colors 
  ), 
  show_annotation_name = TRUE, 
  annotation_name_gp = grid::gpar(fontsize = 14, 
                                  color = "black", 
                                  fontface = "bold"), 
  annotation_label = c("Cytokine", "Pre-Tx"), 
  show_legend = TRUE
) 

#  Row Annotation
my_gene_col <- my_gene_col[rownames(mat), , drop = FALSE] 

row_ha_left <- ComplexHeatmap::rowAnnotation(
  cluster = my_gene_col$Cluster,
  annotation_legend_param = list(
    cluster = list(title = "Cluster", direction = "horizontal",
                   title_gp = grid::gpar(fontsize = 14, color = "black", fontface = "bold"),
                   title_position = "leftcenter", nrow = 1,
                   legend_gp = grid::gpar(fontsize = 14, color = "black", fontface = "plain"))
  ),
  col = list(cluster = row_cluster_colors),
  annotation_name_rot = 0,
  show_annotation_name = FALSE,
  show_legend = TRUE,
  annotation_label = c("Cluster")
)

# get an idea of the ranges of the matrix
q <- quantile(mat, c(0.01, 0.1, 0.5, 0.9, 0.99), na.rm = TRUE)
print(q)

limit <- ceiling(max(abs(mat), na.rm = TRUE)) 

# Define the color function with the specified gradient, focused on 1%/99% quantiles
if(limit < ceiling(abs(min(mat, na.rm = TRUE)))) {
  message("Heatmap limit corresponds to the lower bounds of the logFC range.")
  col_fun = circlize::colorRamp2(c(-limit, 
                                   quantile(x = mat, probs = 0.01, na.rm = TRUE), 
                                   0, 
                                   -quantile(x = mat, probs = 0.01, na.rm = TRUE), 
                                   limit 
  ), 
  c("navy", "navy", "#ededed", "firebrick", "firebrick"))
} else { 
  message("Heatmap limit corresponds to the upper bounds of the logFC range.")
  col_fun = circlize::colorRamp2(c(-limit, 
                                   -quantile(x = mat, probs = 0.99, na.rm = TRUE), 
                                   0, 
                                   quantile(x = mat, probs = 0.99, na.rm = TRUE), 
                                   limit 
  ), 
  c("navy", "navy", "#ededed", "firebrick", "firebrick"))
} 

heatmap_fig <- ComplexHeatmap::Heatmap(
  matrix = mat, 
  heatmap_legend_param=list(
    title="z-score", legend_direction = "horizontal",  
    title_gp = grid::gpar(fontsize = 14, color = "black", fontface = "bold"), 
    labels_gp = grid::gpar(fontsize = 14, color = "black", fontface = "plain") 
  ),  
  na_col = "white", 
  col=col_fun, 
  row_names_gp = grid::gpar(fontsize = 0),
  cluster_rows = hc_row, 
  cluster_columns = hc_col, 
  row_split = k_row, 
  column_title_gp = grid::gpar(fontsize = 0), 
  row_title_gp = grid::gpar(fontsize = 0), 
  top_annotation = column_ha, 
  left_annotation = row_ha_left, 
  show_row_names = FALSE, 
  show_column_names = FALSE, 
  show_row_dend = TRUE, 
  row_dend_side  = "left",
  show_column_dend = FALSE, 
  use_raster = FALSE
)

# Save Heatmap
outdir_heatmap <- file.path(working_directory, "Heatmap")
dir.create(outdir_heatmap, showWarnings = FALSE, recursive = TRUE)

# PDF
grDevices::pdf(file = file.path(outdir_heatmap, paste("differential-RNAseq", "ComplexHeatmap.pdf", sep = "_")),
               height = heatmap.height / 2.54,
               width = heatmap.width / 2.54 
) 

grid::grid.newpage() 
ComplexHeatmap::draw(heatmap_fig, heatmap_legend_side = heatmap.legend.position, annotation_legend_side = heatmap.legend.position, legend_grouping = heatmap.legend.grouping, 
                     padding = unit(c(2, 2, 2, 2), "mm"), use_raster = FALSE) 
dev.off()

# PNG
grDevices::png(filename = file.path(outdir_heatmap, paste("differential-RNAseq", "ComplexHeatmap.png", sep = "_")),
               height = heatmap.height,
               width = heatmap.width,
               res = plot.dpi, 
               units = "cm"
) 

grid::grid.newpage() 
ComplexHeatmap::draw(heatmap_fig, heatmap_legend_side = heatmap.legend.position, annotation_legend_side = heatmap.legend.position, legend_grouping = heatmap.legend.grouping, 
                     padding = unit(c(2, 2, 2, 2), "mm"), use_raster = FALSE) 
dev.off() 

# SVG
grDevices::svg(filename = file.path(outdir_heatmap, paste("differential-RNAseq", "ComplexHeatmap.svg", sep = "_")),
               height = heatmap.height / 2.54,
               width = heatmap.width / 2.54 
) 

grid::grid.newpage() 
ComplexHeatmap::draw(heatmap_fig, heatmap_legend_side = heatmap.legend.position, annotation_legend_side = heatmap.legend.position, legend_grouping = heatmap.legend.grouping, 
                     padding = unit(c(2, 2, 2, 2), "mm"), use_raster = FALSE) 
dev.off() 

print(paste("Number of genes in global DE heatmap", ":", nrow(mat), sep = " "))
 
# ---- Cluster-based transcription factor (TF) enrichment analysis with Discriminant Regulon Expression Analysis (DoRothEA) ---- 

# Convert viper_regulon to TERM2GENE format required by clusterprofiler::enricher()
TERM2GENE <- purrr::map_dfr(
  .x = names(viper_regulon),
  .f = function(tf) {
    targets <- names(viper_regulon[[tf]]$tfmode)
    if(length(targets) > 0) {
      data.frame(
        term = tf,
        gene = targets,
        stringsAsFactors = FALSE
      )
    }
  }
)

# Strip the " - omnipath_scores/TOP/network.sif" suffix from TF names, to make output cleaner
TERM2GENE$term <- sub(" - .*", "", TERM2GENE$term) 

# Keep only tiers A, B, and C (highest confidence)
TERM2GENE <- TERM2GENE %>% dplyr::filter(grepl("_[ABC]$", term))

head(TERM2GENE) 

# Identify background (universe) gene set: 
universe_hgnc <- gene_metadata$gene_symbol %>% na.omit() %>% unique() 
 
# Perform pathway enrichment analysis for each cluster 
dorothea_results <- list() 
for(i in seq_along(cluster.list)) {
  x <- cluster.list[i]
  message (paste0("Starting DoRothEA enrichment analysis for cluster: ", x)) 
  
  # Set output directory to save results
  outdir_tf <- file.path(working_directory, "TF_enrichment_DoRothEA", paste0("cluster", x))
  dir.create(outdir_tf, showWarnings = FALSE, recursive = TRUE)
  
  # Filter genes to elements in the cluster
  query_hgnc <- my_gene_col %>% 
    dplyr::filter(as.integer(as.character(Cluster)) == x) %>% 
    dplyr::pull(gene_symbol) %>% 
    na.omit() %>% unique() 
  
  # Perform TF specific over-representation analysis (ORA) with clusterProfiler::enricher()
  tf_enrich <- clusterProfiler::enricher(
    gene = query_hgnc,
    universe = universe_hgnc, 
    TERM2GENE = TERM2GENE,
    pvalueCutoff = p.value.threshold,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.1,
    minGSSize = 5,
    maxGSSize = 500, 
    gson = NULL
  ) 
  
  if(is.null(tf_enrich)) {
    warning(paste0("clusterProfiler::enricher() returned NULL for cluster ", x, " (no genes could be mapped). Moving to next cluster.")) 
    dorothea_results[[as.character(x)]] <- NA
    next()
  }
  
  if(nrow(as.data.frame(tf_enrich@result)) < 1) {
    warning(paste0("No significantly enriched transcription factors identified in cluster ", x, ". Saving enrichment results and moving to next cluster.")) 
    dorothea_results[[as.character(x)]] <- NA 
    
    # Save the TF enrichment results
    utils::write.csv(x = as.data.frame(tf_enrich@result), file = file.path(outdir_tf, paste("DoRothEA", paste0("cluster", x), "enrichResult.csv", sep = "_")),   row.names = FALSE) 
    
    # Skip the rest of the loop and move to the next iteration (cluster)
    next()
  } 
  
  tf_geneSets <- data.frame(
    ID = names(tf_enrich@geneSets), 
    GeneSet_count = purrr::map_int(tf_enrich@geneSets, length)
  ) 
  
  # Extract results and calculate odds ratio
  tf_enrich_results <- as.data.frame(tf_enrich@result) %>% 
    dplyr::left_join(tf_geneSets, by = "ID") %>% 
    dplyr::mutate(
      universe_count = as.numeric(sub(".*/", "", tf_enrich@result$BgRatio[1])), 
      DEG_count = as.numeric(sub(".*/", "", GeneRatio)), 
      percent.overlap = (Count/GeneSet_count)*100
    ) %>% 
    dplyr::mutate(
      a = Count, 
      b = DEG_count - Count, 
      c = GeneSet_count - Count, 
      d = (universe_count - DEG_count) - (GeneSet_count - Count)
    ) 
  
  # Calculate fisher exact test and perform BH p-value correction
  tf_enrich_results <- tf_enrich_results %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate( 
      fisher = list(stats::fisher.test(matrix(c(a, b, c, d), nrow = 2))), 
      odds_ratio.exact = as.numeric(fisher$estimate), 
      conf.low = fisher$conf.int[1], 
      conf.high = fisher$conf.int[2], 
      p_fisher = fisher$p.value 
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(
      p.adjust_fisher = stats::p.adjust(p_fisher, 
                                        method = "BH")
    ) %>% 
    dplyr::select(-fisher, -a, -b, -c, -d) 
  
  # Separate the DoRothEA confidence level and the TF name into separate columns
  tf_enrich_results <- tf_enrich_results %>% 
    tidyr::separate(Description, into = c("Description", "DoRothEA.confidence"), sep = "_(?=[^_]+$)", remove = TRUE)
  
  # Identify top enriched TFs 
  pathways_up <- tf_enrich_results %>% 
    dplyr::filter(!is.na(p.adjust_fisher) & p.adjust_fisher < p.value.threshold & Count >= count.threshold) %>% 
    dplyr::mutate(primary_rank = -log10(p.adjust_fisher)) %>%   
    dplyr::arrange(dplyr::desc(primary_rank), odds_ratio.exact) %>% 
    dplyr::slice_head(n = 6) %>% 
    dplyr::pull(ID) %>% 
    unique() 
  
  if(length(pathways_up) < 1) {
    warning(paste0("No transcription factors passed the significance + count threshold in cluster ", x, 
                   ". Saving enrichment results and skipping plot.")) 
    dorothea_results[[as.character(x)]] <- NA 
    utils::write.csv(tf_enrich_results, file = file.path(outdir_tf, paste("DoRothEA", paste0("cluster", x), "enrichResult.csv", sep = "_")), row.names = FALSE)
    next()
  }
  
  # Prep dataframe for plotting 
  df <- tf_enrich_results %>%
    dplyr::filter(!is.na(p.adjust_fisher) & ID %in% pathways_up) %>%
    dplyr::mutate(primary_rank = -log10(p.adjust_fisher)) %>%
    dplyr::arrange(primary_rank, dplyr::desc(odds_ratio.exact)) %>%
    dplyr::mutate(
      Description = factor(Description, levels = unique(Description))
    ) 
  
  #Compute relevant scalers for plot
  max_or <- ceiling(max(df$odds_ratio.exact, na.rm = TRUE))
  max_adjP <- ceiling(max(-log10(df$p.adjust_fisher), na.rm = TRUE)) 
  max_overlap <- max(50, ceiling(max(df$percent.overlap, na.rm = TRUE)))
  size_range <- c(3, 12) 
  size_breaks <- c(10, round((max_overlap-10) / 2), max_overlap)
  
  #Draw plot
  plot <- ggplot2::ggplot(data = df) + 
    ggplot2::geom_segment(ggplot2::aes(y = Description, yend = Description, x = 0, xend = -log10(p.adjust_fisher)), color = "black") + 
    ggplot2::geom_point(ggplot2::aes(y = Description, x = -log10(p.adjust_fisher), fill = odds_ratio.exact, size = percent.overlap), shape = 21, color = "black") + 
    ggplot2::scale_fill_gradient(low = "#ededed", high = "darkred", limits = c(1, max_or), name = "Odds Ratio") + 
    ggplot2::scale_x_continuous(limits = c(0, max_adjP), expand = ggplot2::expansion(mult = c(0, 0.15))) + 
    ggplot2::scale_size_continuous(limits = c(1, max_overlap), 
                                   range = size_range, 
                                   breaks = size_breaks, 
                                   labels = paste0(size_breaks, "%"), 
                                   name = "% Overlap") + 
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) + 
    ggplot2::guides(
      fill = ggplot2::guide_colorbar(order = 1),  # color legend first (top)
      size = ggplot2::guide_legend(order = 2)     # size legend second (bottom) 
    ) + 
    ggplot2::labs(x = "-log10(adj. p)", 
                  y = NULL) + 
    ggplot2::theme_minimal() + 
    ggplot2::theme(legend.position = "right", 
                   axis.text.x = ggplot2::element_text(size = 12, 
                                                       face = "plain", 
                                                       color = "black"),  
                   axis.text.y = ggplot2::element_text(size = 14, 
                                                       color = "black", 
                                                       face = "bold"), 
                   axis.title = ggplot2::element_text(size = 14, 
                                                      color = "black", 
                                                      face = "bold"), 
                   legend.text = ggplot2::element_text(size = 12, 
                                                       face = "plain", 
                                                       color = "black"), 
                   legend.title = ggplot2::element_text(size = 14, 
                                                        color = "black", 
                                                        face = "bold")) 
  print(plot) 
  
  # Save the DoRothEA enrichment results
  utils::write.csv(x = tf_enrich_results, file = file.path(outdir_tf, paste("DoRothEA", paste0("cluster", x), "enrichResult.csv", sep = "_")), row.names = FALSE) 
  
  # Save plot
  
  # PDF
  ggplot2::ggsave( 
    paste("ggplot", "DoRothEA", "enrichment", paste0("cluster", x), "LollipopPlot.pdf", sep = "_"), 
    path = outdir_tf, 
    plot = plot, 
    device = cairo_pdf, 
    dpi = plot.dpi, 
    height = lollipop.height, 
    width = lollipop.width, 
    units = "cm", 
    bg = "transparent", 
    limitsize = TRUE) 
  
  # PNG
  ggplot2::ggsave( 
    paste("ggplot", "DoRothEA", "enrichment", paste0("cluster", x), "LollipopPlot.png", sep = "_"), 
    path = outdir_tf, 
    plot = plot, 
    device = "png", 
    dpi = plot.dpi, 
    height = lollipop.height, 
    width = lollipop.width, 
    units = "cm", 
    bg = "transparent", 
    limitsize = TRUE) 
  
  # SVG
  ggplot2::ggsave( 
    paste("ggplot", "DoRothEA", "enrichment", paste0("cluster", x), "LollipopPlot.svg", sep = "_"), 
    path = outdir_tf, 
    plot = plot, 
    device = "svg", 
    dpi = plot.dpi, 
    height = lollipop.height, 
    width = lollipop.width, 
    units = "cm", 
    bg = "transparent", 
    limitsize = TRUE) 
  
  # Save results to DoRothEAR results list
  dorothea_results[[as.character(x)]] <- tf_enrich_results
  
  message (paste0("Finished DoRothEA enrichment analysis for cluster: ", x)) 
} 

# ---- Update heatmap to incorporate DoRothEA TF enrichment results ---- 

# Identify hierarchical clustering + cutree-based cluster assignment 
row_cluster_labels <- my_gene_col$Cluster   # cutree-derived factor, aligned to rownames(mat)

# Build tf_text_list keyed cutree cluster IDs 
tf_text_list <- list()
for (cluster.name in levels(row_cluster_labels)) {
  if (is.null(dorothea_results[[cluster.name]]) ||
      (length(dorothea_results[[cluster.name]]) == 1 && is.na(dorothea_results[[cluster.name]]))) {
    tf_text_list[[cluster.name]] <- " "
    next
  }
  df <- dorothea_results[[cluster.name]]
  tfs <- df %>%
    dplyr::filter(p.adjust_fisher < p.value.threshold & Count >= count.threshold) %>%
    dplyr::arrange(p.adjust_fisher) %>%
    dplyr::slice_min(order_by = p.adjust_fisher, n = 10, na_rm = TRUE) %>%
    dplyr::pull(Description) %>%
    unique()
  tf_text_list[[cluster.name]] <- if (length(tfs) == 0) " " else tfs
}
tf_text_list <- tf_text_list[levels(row_cluster_labels)]

# Build a dataframe with cutree cluster assignments as a column
tf_df <- tibble::enframe(tf_text_list, name = "Cluster", value = "text") %>%
  tidyr::unnest(text) %>%
  dplyr::mutate(
    fontsize   = 11,
    col        = row_cluster_colors[Cluster], 
    fontfamily = "sans",
    fontface   = "bold"
  )

# Draw a preliminary split-only heatmap to learn ComplexHeatmap's internal position order 
ch_row_order <- ComplexHeatmap::row_order(ComplexHeatmap::draw(heatmap_fig))

position_to_cutree <- sapply(seq_along(ch_row_order), function(pos) {
  gene_indices  <- ch_row_order[[pos]]
  gene_names    <- rownames(mat)[gene_indices]
  cutree_labels <- clusters_row[gene_names]
  unique_label  <- unique(cutree_labels)
  if (length(unique_label) != 1) {
    stop(sprintf("Position %d maps to multiple cutree clusters: %s",
                 pos, paste(unique_label, collapse = ", ")))
  }
  as.character(unique_label)
})
names(position_to_cutree) <- as.character(seq_along(position_to_cutree))

# Invert the cluster mapping: cutree ID -> CH position label 
cutree_to_position <- setNames(names(position_to_cutree), position_to_cutree)

# Recode row_cluster_labels into CH's position numberingv
align_to_ch <- factor(
  cutree_to_position[as.character(row_cluster_labels)],
  levels = as.character(seq_along(position_to_cutree))
)

# Build text_list keyed by CH posittion, pulling text/color from tf_df via cutree ID 
text_list <- lapply(names(position_to_cutree), function(pos.name) {
  cutree_id <- position_to_cutree[[pos.name]]
  tf_df %>%
    dplyr::filter(Cluster == cutree_id) %>%
    dplyr::select(text, fontsize, col, fontfamily, fontface)
})
names(text_list) <- names(position_to_cutree)

# Map new row annotation 
row_ha_right <- ComplexHeatmap::rowAnnotation(
  TF = ComplexHeatmap::anno_textbox(
    align_to      = align_to_ch, 
    text          = text_list, 
    max_width     = grid::unit(2, "cm"),
    background_gp = grid::gpar(fill = NA, col = NA) 
  ),
  annotation_legend_param = list(),
  annotation_name_rot = 0,
  show_annotation_name = FALSE,
  show_legend = FALSE
)

# Draw final heatmap ----
heatmap_tf_fig <- ComplexHeatmap::Heatmap( 
  matrix = mat,
  heatmap_legend_param = list(
    title = "z-score", legend_direction = "horizontal",
    title_gp = grid::gpar(fontsize = 14, color = "black", fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 14, color = "black", fontface = "plain")
  ),
  na_col = "white",
  col = col_fun, 
  row_names_gp = grid::gpar(fontsize = 0), 
  cluster_rows = hc_row,
  cluster_columns = hc_col,
  row_split = k_row,  
  column_title_gp = grid::gpar(fontsize = 0), 
  row_title_gp = grid::gpar(fontsize = 0), 
  top_annotation = column_ha,
  left_annotation = row_ha_left,
  right_annotation = row_ha_right,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = TRUE, 
  row_dend_side  = "left", 
  show_column_dend = FALSE,
  use_raster = FALSE
) 

# Save Heatmap
outdir_heatmap <- file.path(working_directory, "Heatmap") 
dir.create(outdir_heatmap, showWarnings = FALSE, recursive = TRUE) 

# PDF
grDevices::pdf(file = file.path(outdir_heatmap, paste("differential-RNAseq", "TF-annotation", "ComplexHeatmap.pdf", sep = "_")),
               height = heatmap.height / 2.54,
               width = (heatmap.width + 2) / 2.54 
) 

grid::grid.newpage() 
ComplexHeatmap::draw(heatmap_tf_fig, heatmap_legend_side = heatmap.legend.position, annotation_legend_side = heatmap.legend.position, legend_grouping = heatmap.legend.grouping, 
                     padding = unit(c(2, 2, 2, 2), "mm"), use_raster = FALSE) 
dev.off()

# PNG
grDevices::png(filename = file.path(outdir_heatmap, paste("differential-RNAseq", "TF-annotation", "ComplexHeatmap.png", sep = "_")),
               height = heatmap.height,
               width = (heatmap.width + 2),
               res = 600, 
               units = "cm"
) 

grid::grid.newpage() 
ComplexHeatmap::draw(heatmap_tf_fig, heatmap_legend_side = heatmap.legend.position, annotation_legend_side = heatmap.legend.position, legend_grouping = heatmap.legend.grouping, 
                     padding = unit(c(2, 2, 2, 2), "mm"), use_raster = FALSE) 
dev.off() 

# SVG
grDevices::svg(filename = file.path(outdir_heatmap, paste("differential-RNAseq", "TF-annotation", "ComplexHeatmap.svg", sep = "_")),
               height = heatmap.height / 2.54,
               width = (heatmap.width + 2) / 2.54
) 

grid::grid.newpage() 
ComplexHeatmap::draw(heatmap_tf_fig, heatmap_legend_side = heatmap.legend.position, annotation_legend_side = heatmap.legend.position, legend_grouping = heatmap.legend.grouping, 
                     padding = unit(c(2, 2, 2, 2), "mm"), use_raster = FALSE) 
dev.off() 

# ---- PROGENy pathway enrichment model per sample ---- 
# Retrieve a log2CPM matrix to input into progeny 
logcpm_matrix <- edgeR::cpm(y_global, log = TRUE) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var = "geneName") %>% 
  dplyr::left_join(gene_metadata, by = "geneName") %>% 
  dplyr::filter(!is.na(gene_symbol), , gene_symbol != "") %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(mean_logcpm = mean(dplyr::c_across(dplyr::starts_with("lib")), na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::slice_max(mean_logcpm, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(gene_symbol, dplyr::starts_with("lib")) %>%
  tibble::column_to_rownames(var = "gene_symbol") %>%
  as.matrix() 

# Checks 
sum(duplicated(rownames(logcpm_matrix)))  # should be 0
anyNA(rownames(logcpm_matrix))            # should be FALSE
any(rownames(logcpm_matrix) == "")        # should be FALSE

# Run pathway enrichment analysis
sample_pathway_scores <- decoupleR::run_mlm(
  mat = logcpm_matrix,
  network = progeny_net,
  .source = "source", .target = "target", .mor = "weight",
  minsize = 5
) %>% 
  dplyr::left_join(sample_metadata %>% dplyr::select(libid, Replicate, cytokine, arm, condition_level), 
                   by = c("condition" = "libid"))

# ---- Descriptive: 3 pathways of interest, contrast-level and sample-level ----
sample_pathway_scores %>%
  dplyr::filter(source %in% c("JAK-STAT", "NFkB", "MAPK")) %>% 
  dplyr::arrange(source, condition) %>%
  print(n = Inf) 

# ---- Formal test: all pathways ---- 
# Build pathway score matrix: pathways (rows) x samples (columns) 
pathway_matrix <- sample_pathway_scores %>%
  tidyr::pivot_wider(id_cols = source, names_from = condition, values_from = score) %>%
  tibble::column_to_rownames("source") %>%
  as.matrix()

# Reorder columns to exactly match design_global's sample order, and check 
rownames(design_global) <- sample_metadata$libid
head(rownames(design_global))

pathway_matrix <- pathway_matrix[, rownames(design_global)]
stopifnot(identical(colnames(pathway_matrix), rownames(design_global)))

# ---- Fit with the blocked design ---- 
# Same JakI - STIM contrast, per cytokine ----
stim_vs_jaki_contrast <- limma::makeContrasts(
  IFNl3 = IFNl3_JakI - IFNl3_Stim,
  IL19  = IL19_JakI - IL19_Stim,
  IL20  = IL20_JakI - IL20_Stim,
  IL22  = IL22_JakI - IL22_Stim,
  IL24  = IL24_JakI - IL24_Stim,
  IL26  = IL26_JakI - IL26_Stim,
  levels = design_global
)

fit_pathway <- limma::lmFit(pathway_matrix, design_global)
fit_pathway_2 <- limma::contrasts.fit(fit_pathway, stim_vs_jaki_contrast)
fit_pathway_2_eb <- limma::eBayes(fit_pathway_2)

# Extract all 6 cytokine results
pathway_stim_vs_jaki <- purrr::map_dfr(
  colnames(stim_vs_jaki_contrast),
  function(cytokine_name) {
    limma::topTable(fit_pathway_2_eb, coef = cytokine_name, number = Inf) %>%
      tibble::rownames_to_column("pathway") %>%
      dplyr::mutate(cytokine = cytokine_name, .before = 1)
  }
) %>% 
  rstatix::add_significance(p.col = "adj.P.Val") %>% 
  dplyr::mutate( 
    cytokine.names = dplyr::case_when(
      cytokine == "IFNl3" ~ "IFNλ3", 
      cytokine == "IL19" ~ "IL-19", 
      cytokine == "IL20" ~ "IL-20", 
      cytokine == "IL22" ~ "IL-22", 
      cytokine == "IL24" ~ "IL-24", 
      cytokine == "IL26" ~ "IL-26", 
      TRUE ~ cytokine
    ), 
    pathway.names = dplyr::case_when(
      pathway == "TNFa" ~ "TNFα", 
      pathway == "TGFb" ~ "TGFβ", 
      pathway == "NFkB" ~ "NF-κB", 
      TRUE ~ pathway
    ), 
    Significance = dplyr::case_when(
      adj.P.Val < 0.05 ~ "Significant", 
      TRUE ~ "Not Significant"
    ) 
  ) %>% 
  dplyr::rename(
    activity.score.difference = logFC, 
    mean.activity.score = AveExpr, 
    log.odds.differential.activity = B
  ) 

head(pathway_stim_vs_jaki)

pathway_stim_vs_jaki %>% 
  dplyr::select(cytokine, pathway, activity.score.difference, t, P.Value, adj.P.Val, Significance) %>%
  dplyr::arrange(adj.P.Val) %>% 
  tibble::as_tibble() %>% 
  print(n = Inf) 

# ---- Draw bar plots of pathway enrichment results (per pathway) ---- 
# Define output directory
outdir_bar <- file.path(working_directory, "PROGENy_pathway_enrichment_bar_plots")
dir.create(outdir_bar, showWarnings = FALSE, recursive = TRUE)

# Define pathways to plot
pathways <- unique(sample_pathway_scores$source) 

for (i in 1:length(pathways)) {
  current.pathway <- pathways[i]
  message(paste0("Plotting ", current.pathway, " pathway."))
  
  # Define mock (background) pathway score
  average.mock.score <- sample_pathway_scores %>%
    dplyr::filter(cytokine == "Mock" & source == current.pathway) %>%
    dplyr::pull(score) %>% mean(na.rm = TRUE)
  
  # Define pathway significance between cytokine STIM +/- panJAKi
  pathway.significance <- pathway_stim_vs_jaki %>%
    dplyr::filter(pathway == current.pathway) %>%
    dplyr::select(cytokine, pathway, Significance, adj.P.Val.signif) %>%  
    dplyr::mutate(
      cytokine = dplyr::if_else(cytokine == "IFNl3", "IFNl3", sub("^IL", "IL-", cytokine)),
      cytokine = factor(cytokine, levels = cytokine_order)
    ) 
  stopifnot(!anyNA(pathway.significance$cytokine))
  
  if (all(pathway.significance$Significance == "Not Significant")) {
    message(paste0(
      "No significant differential regulation for ", current.pathway,
      " between panJAKi pretreated and untreated samples post cytokine stimulation. Skipping to the next pathway."
    ))
    next()
  } else {
    message(paste0(
      "Significant differential regulation is present between panJAKi pretreated and untreated samples post cytokine stimulation for ",
      current.pathway, "."
    ))
  } 
  
  # Generate the plotting data frames 
  plot.df <- sample_pathway_scores %>% 
    dplyr::filter(source == current.pathway & cytokine != "Mock") %>% 
    dplyr::mutate(
      cytokine = factor(cytokine, levels = cytokine_order), 
      arm = factor(arm, levels = c("Stim", "JakI"))
    ) 
  stopifnot(!anyNA(plot.df$arm), !anyNA(plot.df$cytokine))
  
  plot.summary.df <- sample_pathway_scores %>% 
    dplyr::filter(source == current.pathway & cytokine != "Mock") %>% 
    dplyr::group_by(cytokine, arm, condition_level) %>% 
    dplyr::summarise(
      count = dplyr::n(),
      median.score = median(score, na.rm = TRUE),
      upper.quartile = stats::quantile(score, probs = 0.75, na.rm = TRUE), 
      lower.quartile = stats::quantile(score, probs = 0.25, na.rm = TRUE),
      .groups = "drop"
    ) %>% 
    dplyr::mutate(
      cytokine = factor(cytokine, levels = cytokine_order), 
      arm = factor(arm, levels = c("Stim", "JakI")) 
    ) 
  stopifnot(!anyNA(plot.summary.df$arm), !anyNA(plot.summary.df$cytokine)) 
  
  # Calculate significance label height, then attach to pathway.significance 
  data_max <- max(c(plot.df$score, plot.summary.df$upper.quartile), na.rm = TRUE)
  data_min <- min(c(plot.df$score, plot.summary.df$lower.quartile), na.rm = TRUE)
  data_range <- data_max - data_min
  visual_top <- max(0, data_max) 
  upper_limit <- visual_top + data_range * 0.4   # reserve 40% of the data range as headroom, once
  segment_y <- visual_top + data_range * 0.14    # sits low in the reserved margin
  
  pathway.significance <- pathway.significance %>%
    dplyr::mutate(segment_y = segment_y)
  
  # Generate cytokine labels 
  cytokine_labels <- pathway_stim_vs_jaki %>% 
    dplyr::mutate(
    cytokine = dplyr::if_else(cytokine == "IFNl3", "IFNl3", sub("^IL", "IL-", cytokine))
  ) %>%
    dplyr::distinct(cytokine, cytokine.names) %>% 
    dplyr::arrange(cytokine) %>% 
    tibble::deframe() 
  
  # Generate faceted bar plots 
  bar.plot <- ggplot2::ggplot(data = NULL, ggplot2::aes(x = arm)) + 
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = "black") + 
    ggplot2::geom_hline(yintercept = average.mock.score, linewidth = 0.5, linetype = "dashed", color = "grey40") + 
    ggplot2::geom_errorbar(data = plot.summary.df, ggplot2::aes(ymin = median.score, ymax = upper.quartile), color = "black", linewidth = 0.4) + 
    ggplot2::geom_errorbar(data = plot.summary.df, ggplot2::aes(ymin = lower.quartile, ymax = median.score), color = "black", linewidth = 0.4) + 
    ggplot2::geom_bar(data = plot.summary.df, stat = "identity", ggplot2::aes(y = median.score, color = arm), linewidth = 0.7, fill = rgb(1, 0, 0, alpha = 0)) + 
    ggplot2::geom_line(
      data = plot.df, 
      ggplot2::aes(y = score, group = Replicate), 
      color = "grey40", linewidth = 0.4, alpha = 0.9
    ) + 
    ggplot2::geom_point(
      data = plot.df,
      ggplot2::aes(y = score, fill = arm),
      shape = 21, color = "black", size = 2.5, alpha = 0.9
    ) + 
    ggplot2::geom_segment(
      data = pathway.significance,
      ggplot2::aes(x = 1, xend = 2, y = segment_y, yend = segment_y),
      inherit.aes = FALSE, linewidth = 0.5, color = "black"
    ) + 
    ggplot2::geom_text(
      data = pathway.significance,
      ggplot2::aes(x = 1.5, y = Inf, label = adj.P.Val.signif),
      inherit.aes = FALSE, size = 4, fontface = "bold", color = "black"
    ) + 
    ggplot2::scale_fill_manual(values = c(Stim = "black", "JakI" = "magenta3"), 
                               labels = c(Stim = "Mock", "JakI" = "panJAKi")) + 
    ggplot2::scale_color_manual(values = c(Stim = "black", "JakI" = "magenta3"), 
                                labels = c(Stim = "Mock", "JakI" = "panJAKi")) + 
    ggplot2::theme_minimal() + 
    ggplot2::theme(
      panel.spacing = ggplot2::unit(0, "lines"), 
      axis.text.y = ggplot2::element_text(color = "black", size = 9, face = "plain"), 
      axis.text.x = ggplot2::element_blank(), 
      axis.title.y = ggplot2::element_text(color = "black", size = 10, face = "bold"), 
      axis.title.x = ggplot2::element_blank(), 
      legend.text = ggplot2::element_text(color = "black", size = 9, face = "plain"), 
      legend.title = ggplot2::element_text(color = "black", size = 10, face = "bold"), 
      strip.text = ggplot2::element_text(color = "black", size = 10, face = "bold", angle = 45), 
      legend.key.size = ggplot2::unit(0.4, "cm"),      # smaller legend keys/swatches
      legend.spacing.y = ggplot2::unit(0.2, "cm"),     # tighter vertical spacing between legends
      legend.margin = ggplot2::margin(0, 0, 0, 0),     # remove legend internal padding
      legend.box.spacing = ggplot2::unit(0.2, "cm")    # tighter gap between plot and legend
      
    ) + 
    ggplot2::coord_cartesian(clip = "off") + 
    ggplot2::labs(title = paste0("Pathway: ", current.pathway), y = "PROGENy\nPathway Activity Score", 
                  color = "Pretreatment", fill = "Pretreatment") + 
    ggplot2::facet_wrap(~ cytokine, nrow = 1, scales = "fixed", strip.position = "bottom", 
                        labeller = ggplot2::as_labeller(cytokine_labels)) 
  
  # Save plots 
  
  # PDF
  ggplot2::ggsave(
    filename = paste("RNAseq", "JAKi", "PROGENy", current.pathway, "PathwayEnrichment", "barplot", "pdf", sep = "."), 
    path = outdir_bar, 
    plot = bar.plot, 
    device = cairo_pdf, 
    dpi = plot.dpi, 
    height = bar.height, 
    width = bar.width, 
    units = "cm", 
    limitsize = TRUE) 
  
  # PNG 
  ggplot2::ggsave(
    filename = paste("RNAseq", "JAKi", "PROGENy", current.pathway, "PathwayEnrichment", "barplot", "png", sep = "."), 
    path = outdir_bar, 
    plot = bar.plot, 
    device = "png", 
    dpi = plot.dpi, 
    height = bar.height, 
    width = bar.width, 
    units = "cm", 
    limitsize = TRUE) 
  
  # SVG 
  ggplot2::ggsave(
    filename = paste("RNAseq", "JAKi", "PROGENy", current.pathway, "PathwayEnrichment", "barplot", "svg", sep = "."), 
    path = outdir_bar, 
    plot = bar.plot, 
    device = "svg", 
    dpi = plot.dpi, 
    height = bar.height, 
    width = bar.width, 
    units = "cm", 
    limitsize = TRUE) 
  
  message(paste0("Finishing ", current.pathway, " pathway."))
}  

# ---- Dot plot summary plot of pathways of interest (JAK-STAT, MAPK, NFkB) 
# Define pathways of interest 
pathways.of.interest <- c("JAK-STAT", "MAPK", "NFkB") 

# Prep the data object 
dot.plot.df <- pathway_stim_vs_jaki %>% 
  dplyr::filter(pathway %in% pathways.of.interest) %>% 
  dplyr::mutate( 
    cytokine = dplyr::if_else(cytokine == "IFNl3", "IFNl3", sub("^IL", "IL-", cytokine)), 
    cytokine = factor(cytokine, levels = cytokine_order), 
    pathway = factor(pathway, levels = pathways.of.interest)
  ) 
stopifnot(!anyNA(dot.plot.df$cytokine), !anyNA(dot.plot.df$pathway)) 

# Specify plotting parameters 
delta.max <- max(abs(dot.plot.df$activity.score.difference), na.rm = TRUE)
med <- median(dot.plot.df$activity.score.difference, na.rm = TRUE)
mad_val <- stats::mad(dot.plot.df$activity.score.difference, na.rm = TRUE, constant = 1.4826)

# Compress to specified median absolute deviation for visualization. 
compress_at <- min(mad.limit * mad_val, delta.max * 0.5, na.rm = TRUE)  # safety cap so it never exceeds half the range

breaks <- c(-delta.max, -compress_at, 0, compress_at, delta.max)

# scale_fill_gradientn's `values` must be rescaled to [0, 1], relative to `limits`
rescaled_values <- scales::rescale(breaks, to = c(0, 1), from = c(-delta.max, delta.max)) 

delta.p.limit <- ceiling(max(-log10(dot.plot.df$adj.P.Val), na.rm = TRUE) * 10) / 10 

# Define labels 
row_labels <- dot.plot.df %>% 
  dplyr::distinct(pathway, pathway.names) %>% 
  dplyr::arrange(pathway) %>%
  tibble::deframe() 

cytokine_labels <- dot.plot.df %>% 
  dplyr::distinct(cytokine, cytokine.names) %>% 
  dplyr::arrange(cytokine) %>% 
  tibble::deframe() 

# Mark contrasts that are visually capped (fill) to add annotation
dot.plot.df <- dot.plot.df %>% 
  dplyr::mutate(
    capped = dplyr::case_when(
      abs(activity.score.difference) > compress_at ~ TRUE, 
      TRUE ~ FALSE
    )
  )

# Draw the dotplot 
pathway.dot.plot <- ggplot2::ggplot(data = dot.plot.df, 
                                    ggplot2::aes(x = cytokine, y = pathway, fill = activity.score.difference, size = -log10(adj.P.Val), color = Significance)
) + 
  ggplot2::geom_point(stat = "identity", shape = 21, stroke = 1) + 
  ggplot2::geom_text(data = dplyr::filter(dot.plot.df, capped == TRUE), 
                     ggplot2::aes(x = cytokine, y = pathway, label = "\u2020"), 
                     color = "black", size = 4, fontface = "bold", vjust = 0.2, inherit.aes = FALSE  
  ) + 
  ggplot2::scale_x_discrete(labels = cytokine_labels) + 
  ggplot2::scale_y_discrete(labels = row_labels, limits = rev) + 
  ggplot2::scale_fill_gradientn(
    colors = c("lightblue", "lightblue", "black", "yellow2", "yellow2"),
    values = rescaled_values,
    limits = c(-delta.max, delta.max)
  ) + 
  ggplot2::scale_color_manual(values = c("Significant" = "magenta3", "Not Significant" = "black")) + 
  ggplot2::scale_size(
    range = c(2, 14),
    transform = "sqrt",
    breaks = c(1.3, 5, 10),   # 1.3 ≈ -log10(0.05), a meaningful reference point
    limits = c(0, ceiling(max(-log10(dot.plot.df$adj.P.Val), na.rm = TRUE)))
  ) + 
  ggplot2::theme_minimal() + 
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black", size = 10, face = "bold"), 
    axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1), 
    axis.title = ggplot2::element_blank(), 
    legend.text = ggplot2::element_text(color = "black", size = 9, face = "plain"), 
    legend.title = ggplot2::element_text(color = "black", size = 10, face = "bold"), 
    legend.key.size = ggplot2::unit(0.4, "cm"),      # smaller legend keys/swatches
    legend.spacing.y = ggplot2::unit(0.2, "cm"),     # tighter vertical spacing between legends
    legend.margin = ggplot2::margin(0, 0, 0, 0),     # remove legend internal padding
    legend.box.spacing = ggplot2::unit(0.2, "cm")    # tighter gap between plot and legend
  ) + 
  ggplot2::labs(
    size = "-log10(p value)", 
    color = "Significance", 
    fill = "Δ PROGENy pathway activity\nscore (panJAKi - Mock)"
  ) + 
  ggplot2::guides(
    fill = ggplot2::guide_colorbar(
      order = 1, 
      direction = "horizontal", 
      title.position = "top", 
      barwidth = ggplot2::unit(4.5, "cm"), 
      barheight = ggplot2::unit(0.4, "cm"), 
      angle = 45
    ), 
    color = ggplot2::guide_legend(
      order = 2,
      override.aes = list(size = 3) 
    ), 
    size = ggplot2::guide_legend(
      order = 3, 
      direction = "horizontal", 
      title.position = "top", 
      label.position = "bottom"
    ) 
  ) 

print(pathway.dot.plot)

# ---- Save dot plot ----  
# Define output directory
outdir_dot <- file.path(working_directory, "PROGENy_pathway_enrichment_dot_plot")
dir.create(outdir_dot, showWarnings = FALSE, recursive = TRUE)

# PDF
ggplot2::ggsave(
  filename = paste("RNAseq", "JAKi", "PROGENy", "PathwayEnrichment", "dotplot", "pdf", sep = "."), 
  path = outdir_dot, 
  plot = pathway.dot.plot, 
  device = cairo_pdf, 
  dpi = plot.dpi, 
  height = dot.height, 
  width = dot.width, 
  units = "cm", 
  limitsize = TRUE) 

# PNG 
ggplot2::ggsave(
  filename = paste("RNAseq", "JAKi", "PROGENy", "PathwayEnrichment", "dotplot", "png", sep = "."), 
  path = outdir_dot, 
  plot = pathway.dot.plot, 
  device = "png", 
  dpi = plot.dpi, 
  height = dot.height, 
  width = dot.width, 
  units = "cm", 
  limitsize = TRUE) 

# SVG 
ggplot2::ggsave(
  filename = paste("RNAseq", "JAKi", "PROGENy", "PathwayEnrichment", "dotplot", "svg", sep = "."), 
  path = outdir_dot, 
  plot = pathway.dot.plot, 
  device = "svg", 
  dpi = plot.dpi, 
  height = dot.height, 
  width = dot.width, 
  units = "cm", 
  limitsize = TRUE) 

# ---- Save PROGENy enrichment results ---- 
outdir_prog <- file.path(working_directory, "PROGENy_pathway_enrichment_results")
dir.create(outdir_prog, showWarnings = FALSE, recursive = TRUE) 

utils::write.csv(x = pathway_stim_vs_jaki, file = file.path(outdir_prog, "RNAseq_PROGENy_DifferentialActivity_results.csv"), row.names = FALSE) 
utils::write.csv(x = sample_pathway_scores, file = file.path(outdir_prog, "RNAseq_PROGENy_perSample_activityScores.csv"), row.names = FALSE) 
utils::write.csv(x = gene_metadata, file = file.path(outdir_prog, "RNAseq_panJakI_geneAnnotations.csv"), row.names = FALSE)
utils::write.csv(x = sample_metadata, file = file.path(outdir_prog, "RNAseq_panJakI_sampleMetadata.csv"), row.names = FALSE)

# ---- Export session info ---- 
writeLines(
  utils::capture.output(utils::sessionInfo()),
  con = file.path(working_directory, "session_info.txt")
) 

print("[SUCCESS] JAKi RNAseq analysis pipeline complete.") 