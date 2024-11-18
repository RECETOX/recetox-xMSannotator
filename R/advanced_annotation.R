create_row <- function(adduct, query) {
  query$Adduct <- adduct
  return(query)
}

append_adduct <- function(..., adduct_names) {
  query <- tibble(...)
  rows_with_adduct <- lapply(
    adduct_names,
    create_row,
    query = query
  )
  return(rows_with_adduct)
}

create_chemCompMZ <- function(database, adduct_names) {
  database <- purrr::pmap_dfr(
    database,
    ~ append_adduct(...,
      adduct_names = adduct_names
    )
  )
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
    filter(!duplicated(.data$chemical_ID)) %>%
    count(.data$Confidence)

  confidence_distribution_across_formulas <- annotation %>%
    filter(!duplicated(.data$Formula)) %>%
    count(.data$Confidence)

  print("Confidence level distribution for unique compounds")
  print(confidence_distribution_across_compounds)
  print("Confidence level distribution for unique formulas")
  print(confidence_distribution_across_formulas)

  invisible(annotation)
}

#' Wrapper for the advanced annotation steps.
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
advanced_annotation <- function(peak_table,
                                compound_table,
                                adduct_table = NULL,
                                adduct_weights = NULL,
                                intensity_deviation_tolerance = 0.1,
                                mass_tolerance = 5e-6,
                                mass_defect_tolerance = 0.1,
                                mass_defect_precision = 0.01,
                                time_tolerance = 10,
                                peak_rt_width = 1,
                                correlation_threshold = 0.7,
                                deep_split = 2,
                                min_cluster_size = 10,
                                maximum_isotopes = 10,
                                min_ions_per_chemical = 2,
                                filter_by = c("M-H", "M+H"),
                                network_type = "unsigned",
                                redundancy_filtering = TRUE,
                                outloc = tempdir(),
                                n_workers = parallel::detectCores()) {
  if (is.null(adduct_table)) {
    adduct_table <- sample_adduct_table
  }

  if (is.null(adduct_weights)) {
    adduct_weights <- as.data.frame(tibble(adduct = adduct_table$adduct, weight = rep_len(5, length(adduct_table$adduct))))
  }

  if (is.numeric(n_workers) && n_workers > 1) {
    WGCNA::allowWGCNAThreads(n_workers)
  }

  peak_table <- as_peak_table(peak_table, intensities = TRUE)
  adduct_table <- as_adduct_table(adduct_table)
  compound_table <- as_compound_table(compound_table)

  # Tool 1: Simple annotation
  # ---------------------------
  annotation <- simple_annotation(
    peak_table = peak_table,
    compound_table = compound_table,
    adduct_table = adduct_table,
    mass_tolerance = mass_tolerance
  )

  # this filter is already done in simple_annotation - maybe it can be removed here?
  annotation <- filter(annotation, forms_valid_adduct_pair(.data$molecular_formula, .data$adduct))
  # ---------------------------
 
  # Tool 2: Compute mass defect
  # ----------------------------
  annotation <- compute_mass_defect(annotation, precision = mass_defect_precision)
  # ----------------------------


  # Tool 3: Compute correlations (peak intensity matrix will be the input, so doesn't need to be cut)
  # also check the galaxy tool to compute correlations
  # ----------------------------
  peak_intensity_matrix <- get_peak_intensity_matrix(peak_table)
  peak_correlation_matrix <- compute_peak_correlations(peak_intensity_matrix, correlation_method = "p")
  # ----------------------------

  # Tool 4: Compute peak modules
  # ----------------------------
  peak_modules <- compute_peak_modules(
    peak_intensity_matrix = peak_intensity_matrix,
    peak_correlation_matrix = peak_correlation_matrix,
    correlation_threshold = correlation_threshold,
    deep_split = deep_split,
    min_cluster_size = min_cluster_size,
    network_type = network_type
  )
  # ----------------------------

  # Tool 5: Compute rt modules
  # ----------------------------
  peak_rt_clusters <- compute_rt_modules(
    peak_table = inner_join(peak_table, peak_modules, by = "peak"),
    peak_width = peak_rt_width
  )
  
  # Step can be done in Galaxy using the cut and join tools
  # ----------------------------
  annotation <- inner_join(annotation,
    select(peak_rt_clusters, "peak", "mean_intensity", "module", "rt_cluster"),
    by = "peak"
  )
  # ----------------------------

  # Re-use the mass defect tool on the output of compute rt modules
  # ----------------------------
  peak_table <- peak_table %>%
    select(peak, mz, rt) %>%
    inner_join(peak_rt_clusters, by = "peak") %>%
    compute_mass_defect(precision = mass_defect_precision)
  # ----------------------------


  # Tool 6: Compute isotopes
  # ----------------------------
  annotation <- compute_isotopes(
    annotation = annotation,
    adduct_weights = adduct_weights,
    intensity_deviation_tolerance = intensity_deviation_tolerance,
    mass_defect_tolerance = mass_defect_tolerance,
    peak_table = peak_table,
    rt_tolerance = time_tolerance
  )
  # ----------------------------

  # Tool 7: Reformat annotation and correlation matrix for old xmsannotator
  # ----------------------------
  annotation <- reformat_annotation_table(annotation)
  global_cor <- reformat_correlation_matrix(peak_table, peak_correlation_matrix)
  # ----------------------------

  # Tool 8: Compute chemscores
  # ----------------------------
  annotation <- purrr::pmap_dfr(
    annotation,
    ~ get_chemscore(...,
                    annotation = annotation,
                    adduct_weights = adduct_weights,
                    corthresh = correlation_threshold,
                    global_cor = global_cor,
                    max_diff_rt = time_tolerance,
                    filter.by = filter_by,
                    outlocorig = outloc
    )
  )
  # ----------------------------


  # Tool 9: pathway matching
  # ----------------------------
  data(hmdbCompMZ)
  chemCompMZ <- dplyr::rename(hmdbCompMZ, chemical_ID = HMDBID)

  annotation <- multilevelannotationstep3(
    chemCompMZ = chemCompMZ,
    chemscoremat = annotation,
    adduct_weights = adduct_weights,
    db_name = "HMDB",
    max_diff_rt = time_tolerance,
    pathwaycheckmode = "pm"
  )
  # ----------------------------

  # Tool 10: compute confidence levels
  # ----------------------------
  annotation <- multilevelannotationstep4(
    outloc = outloc,
    chemscoremat = annotation,
    max.mz.diff = mass_tolerance,
    max.rt.diff = time_tolerance,
    filter.by = filter_by,
    adduct_weights = adduct_weights,
    max_isp = maximum_isotopes,
    min_ions_perchem = min_ions_per_chemical
  )
  # ----------------------------

  # Tool 11: print confidence distribution
  # ----------------------------
  print_confidence_distribution(annotation)
  # ----------------------------

  if (redundancy_filtering) {
    # Tool 12: redundancy filtering
    # ----------------------------
    annotation <- multilevelannotationstep5(
      outloc = outloc,
      adduct_weights = adduct_weights,
      chemscoremat = annotation
    )
    # ----------------------------
    
    # Re-use the confidence distrivution printing tool
    # ----------------------------
    print_confidence_distribution(annotation)
    # ----------------------------
  }

  annotation
}
