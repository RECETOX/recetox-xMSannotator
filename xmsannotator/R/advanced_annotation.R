create_chemCompMZ <- function(database, adduct_names) {
  database <- database[which(database$Adduct %in% adduct_names), ]
  return(database)
}

#' @import dplyr
#' @importFrom rlang .data
compute_mass_defect <- function(peaks, precision) {
  mutate(peaks, mass_defect = cut(.data$mz %% 1, seq(0, 1, precision), labels = FALSE))
}

#' @import dplyr
#' @importFrom rlang .data
remove_duplicates <- function(annotation, adduct_weights) {
  is_max <- function(x) x == max(x)
  is_unique <- function(score, adduct) is_max(score * 100^is.element(adduct, adduct_weights$adduct))
  recompute_score <- function(mask, score) score * (sum(mask) / length(mask)) # FIXME: check if decreasing the score is OK

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
#' @importFrom magrittr %>%
advanced_annotation <- function(peak_table,
                                compound_table,
                                pathway_mapping = NULL,
                                exluded_pathways = NULL,
                                exluded_pathway_compounds = NULL,
                                adduct_table = NULL,
                                adduct_weights = NULL,
                                intensity_deviation_tolerance = 0.1,
                                mass_tolerance = 5e-6,
                                mass_defect_tolerance = 0.1,
                                time_tolerance = 10,
                                peak_rt_width = 1,
                                correlation_threshold = 0.7,
                                deep_split = 2,
                                min_cluster_size = 10,
                                network_type = "unsigned",
                                expected_adducts = NULL,
                                boosted_compounds = NULL,
                                redundancy_filtering = TRUE,
                                n_workers = parallel::detectCores()) {
  
  # TODO
  # - these arguments are probably needed: 
  #   db_name, filter.by, pathwaycheckmode (?), max_isp, min_ions_perchem
  # - chemCompMZ needs to be created
  database <- load_database(db_name)
  chemCompMZ <- create_chemCompMZ(database, adduct_names)
  # the adduct_names are constucted from adduct_table and queryadductlist (see multiannot. L139-150)
  # - outloc is still needed in some functions
  outloc = "."
  
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

  annotation <- simple_annotation(
    peak_table = peak_table,
    compound_table = compound_table,
    adduct_table = adduct_table,
    mass_tolerance = mass_tolerance
  )

  peak_modules <- compute_peak_modules(
    peak_intensity_matrix = peak_intensity_matrix,
    peak_correlation_matrix = peak_correlation_matrix,
    correlation_threshold = correlation_threshold,
    deep_split = deep_split,
    min_cluster_size = min_cluster_size,
    network_type = network_type
  )

  peak_rt_clusters <- compute_rt_modules(
    peak_table = inner_join(peak_table, peak_modules, by = "peak"),
    peak_width = peak_rt_width
  )

  peak_table <- peak_table %>%
    select(peak, mz, rt) %>%
    inner_join(peak_rt_clusters, on = "peak") %>%
    compute_mass_defect(precision = 0.01)

  annotation <- filter(annotation, forms_valid_adduct_pair(.data$molecular_formula, .data$adduct))
  annotation <- compute_mass_defect(annotation, precision = 0.01)
  annotation <- inner_join(annotation,
    select(peak_rt_clusters, "peak", "mean_intensity", "module", "rt_cluster"),
    by = "peak"
  )

  annotation <- compute_isotopes(
    annotation = annotation,
    adduct_weights = adduct_weights,
    intensity_deviation_tolerance = intensity_deviation_tolerance,
    mass_defect_tolerance = mass_defect_tolerance,
    peak_table = peak_table,
    rt_tolerance = time_tolerance
  )
  
  annotation <- purrr::pmap_dfr(
    annotation,
    ~ get_chemscore(...,
                    annotation = annotation,
                    adduct_weights = adduct_weights,
                    correlation_threshold = 0.7,
                    peak_correlation_matrix = peak_correlation_matrix,
                    time_tolerance = time_tolerance,
                    outlocorig = outloc
    )
  )

  annotation <- multilevelannotationstep3(
    chemCompMZ = chemCompMZ,
    chemscoremat = annotation,
    adduct_weights = adduct_weights,
    db_name = db_name,
    max_diff_rt = time_tolerance,
    pathwaycheckmode = "pm"
  )
  
  annotation <- multilevelannotationstep4(
    outloc = outloc,
    chemscoremat = annotation,
    max.mz.diff = mass_tolerance,
    max.rt.diff = time_tolerance,
    filter.by = filter.by,
    adduct_weights = adduct_weights,
    max_isp = max_isp,
    min_ions_perchem = min_ions_perchem
  )
  
  annotation <- multilevelannotationstep5(
    outloc = outloc,
    adduct_weights = adduct_weights,
    chemscoremat = annotation
  )

  annotation
}
