# Only place the two designs differ. P560-2 has an inhibitor-only arm, so
# contrasts 3 and 5 exist; P560-4 does not.

expected_conditions <- function(dataset) {
  cfg <- DATASETS[[dataset]]
  arms <- as.vector(t(outer(CYTOKINES, c(cfg$vehicle, cfg$inhibitor), paste, sep = "_")))
  c(cfg$vehicle, if (!is.na(cfg$base_arm)) cfg$base_arm, arms)
}

make_contrasts <- function(dataset) {
  cfg <- DATASETS[[dataset]]
  inh <- cfg$inhibitor; veh <- cfg$vehicle; base <- cfg$base_arm

  dplyr::bind_rows(lapply(CYTOKINES, function(cyt) {
    out <- tibble::tribble(
      ~n, ~expr, ~name, ~label,
      1, sprintf("label%s_%s - label%s", cyt, veh, veh),
         sprintf("%s_1_response_inhibitor_absent", cyt),
         sprintf("%s: cytokine response, %s absent", cyt, inh),
      2, sprintf("label%s_%s - label%s", cyt, inh, veh),
         sprintf("%s_2_net_inhibitor_present", cyt),
         sprintf("%s: net state, cytokine + %s (vs %s)", cyt, inh, veh),
      4, sprintf("label%s_%s - label%s_%s", cyt, inh, cyt, veh),
         sprintf("%s_4_inhibitor_effect_in_stimulated", cyt),
         sprintf("%s: effect of %s, cytokine held constant", cyt, inh)
    )

    if (!is.na(base)) {
      out <- dplyr::bind_rows(out, tibble::tribble(
        ~n, ~expr, ~name, ~label,
        3, sprintf("label%s_%s - label%s", cyt, inh, base),
           sprintf("%s_3_response_within_blocked", cyt),
           sprintf("%s: cytokine response within %s-treated cells", cyt, inh),
        5, sprintf("(label%s_%s - label%s) - (label%s_%s - label%s)",
                   cyt, inh, base, cyt, veh, veh),
           sprintf("%s_5_interaction", cyt),
           sprintf("%s: %s-dependent component of the response", cyt, inh)
      ))
    }
    dplyr::mutate(dplyr::arrange(out, n), cytokine = cyt)
  }))
}

# two columns per cytokine; feeds the all-contrast heatmap
heatmap_contrast_defs <- function(dataset) {
  cfg <- DATASETS[[dataset]]
  dplyr::bind_rows(lapply(CYTOKINES, function(cyt) {
    tibble::tibble(
      cytokine = cyt,
      inhibitor_state = c("absent", "present"),
      expr = c(sprintf("label%s_%s - label%s", cyt, cfg$vehicle, cfg$vehicle),
               sprintf("label%s_%s - label%s", cyt, cfg$inhibitor, cfg$vehicle)))
  })) |>
    dplyr::mutate(contrast = paste(cytokine, inhibitor_state, sep = "_"))
}

require_interaction <- function(dataset) {
  if (!isTRUE(DATASETS[[dataset]]$interaction))
    stop(dataset, " has no inhibitor-only arm; interaction not identifiable.")
  invisible(TRUE)
}
