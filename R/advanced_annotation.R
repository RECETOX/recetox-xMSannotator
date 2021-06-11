#' @import dplyr
#' @importFrom rlang .data
compute_mass_defect <- function (peaks, precision) {
  mutate(peaks, mass_defect = cut(.data$mz %% 1, seq(0, 1, precision), labels = FALSE))
}

#' @import dplyr
#' @importFrom rlang .data
remove_duplicates <- function(annotation, adduct_weights) {
  is_max <- function (x) x == max(x)
  is_unique <- function (score, adduct) is_max(score * 100^is.element(adduct, adduct_weights$adduct))
  recompute_score <- function (mask, score) score * (sum(mask) / length(mask)) # FIXME: check if decreasing the score is OK

  annotation <- group_by(annotation, .data$mz)
  annotation <- mutate(annotation, unique = is_unique(.data$score, .data$adduct))
  annotation <- group_by(annotation, .data$compound)
  annotation <- mutate(annotation, score = recompute_score(.data$unique, .data$score))
  annotation <- ungroup(annotation)
  annotation <- filter(annotation, .data$unique)
  annotation <- select(annotation, -.data$unique, -.data$multiple_match)
  annotation
}

#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
print_confidence_distribution <- function(annotation) {
  confidence_distribution_across_compounds <- annotation %>%
    filter(!duplicated(.data$compound)) %>%
    count(.data$confidence)

  confidence_distribution_across_formulas <- annotation %>%
    filter(!duplicated(.data$molecular_formula)) %>%
    count(.data$confidence)

  print("Confidence level distribution for unique compounds")
  print(confidence_distribution_across_compounds)
  print("Confidence level distribution for unique formulas")
  print(confidence_distribution_across_formulas)

  invisible(annotation)
}

#' @export
#' @import dplyr
advanced_annotation <- function(
  peak_table,
  compound_table,
  pathway_mapping = NULL,
  exluded_pathways = NULL,
  exluded_pathway_compounds = NULL,
  adduct_table = NULL,
  adduct_weights = NULL,
  mass_tolerance = 5e-6,
  time_tolerance = 10,
  correlation_threshold = 0.7,
  deep_split = 2,
  min_cluster_size = 10,
  network_type = "unsigned",
  expected_adducts = NULL,
  boosted_compounds = NULL,
  redundancy_filtering = TRUE,
  n_workers = parallel::detectCores()
) {
  if (is.null(adduct_table)) {
    adduct_table <- sample_adduct_table
  }

  if (is.null(adduct_weights)) {
    adduct_weights <- tibble(adduct = c("M+H", "M-H"), weight = c(5, 5))
  }

  if (is.numeric(n_workers) && n_workers > 1) {
    WGCNA::allowWGCNAThreads(n_workers)
  }

  peak_table <- as_peak_table(peak_table, intensities = TRUE)
  adduct_table <- as_adduct_table(adduct_table)
  compound_table <- as_compound_table(compound_table)

  peak_intensity_matrix <- t(select(peak_table, starts_with("intensity")))
  peak_intensity_matrix <- magrittr::set_colnames(peak_intensity_matrix, peak_table$peak)
  peak_correlation_matrix <- WGCNA::cor(peak_intensity_matrix, use = "p", method = "p")

  peak_modules <- compute_peak_modules(
    peak_intensity_matrix = peak_intensity_matrix,
    peak_correlation_matrix = peak_correlation_matrix,
    correlation_threshold = correlation_threshold,
    deep_split = deep_split,
    min_cluster_size = min_cluster_size,
    network_type = network_type
  )

  annotation <- simple_annotation(
    peak_table = peak_table,
    compound_table = compound_table,
    adduct_table = adduct_table,
    mass_tolerance = mass_tolerance
  )

  annotation <- compute_mass_defect(annotation, precision = 0.01)
  annotation <- inner_join(annotation, peak_modules, by = "peak")

  annotation <- compute_scores(
    annotation = annotation,
    adduct_weights = adduct_weights,
    time_tolerance = time_tolerance
  )

  annotation <- compute_pathways(
    annotation = annotation,
    pathway_mapping = pathway_mapping,
    exluded_pathways = exluded_pathways,
    exluded_pathway_compounds = exluded_pathway_compounds,
    adduct_weights = adduct_weights,
    score_threshold = 0.1
  )

  annotation <- compute_confidence_levels(
    annotation = annotation,
    expected_adducts = expected_adducts,
    boosted_compounds = boosted_compounds,
    mass_tolerance = mass_tolerance,
    time_tolerance = time_tolerance
  )

  print_confidence_distribution(annotation)

  if (redundancy_filtering) {
    annotation <- remove_duplicates(annotation, adduct_weights)
    print_confidence_distribution(annotation)
  }

  annotation
}
