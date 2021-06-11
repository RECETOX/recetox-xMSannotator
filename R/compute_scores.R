#' @import dplyr
#' @importFrom rlang .data
compute_isotopes <- function(peaks, query, rt_tolerance, mass_defect_tolerance = 0, abudance_ratio, max_isp) {
  isotopes <- mutate(
    peaks,
    compound = query$compound,
    isotope_number = round(mz - query$monoisotopic_mass),
    adduct = sprintf("%s_[%+d]", query$adduct, isotope_number),
    formula = sprintf("%s_[%+d]", query$formula, isotope_number),
  )
  isotopes <- filter(
    isotopes,
    cluster == query$cluster,
    avg_intensity / query$avg_intensity <= abudance_ratio,
    near(rt, query$rt, rt_tolerance),
    near(mass_defect, query$mass_defect, mass_defect_tolerance),
    between(abs(isotope_number), 1, max_isp)
  )
  isotopes <- select(isotopes, -isotope_number)
}

#' @import dplyr
#' @importFrom rlang .data
compute_scores <- function(annotation, adduct_weights, time_tolerance) {
  annotation <- filter(annotation, forms_valid_adduct_pair(.data$molecular_formula, .data$adduct))
  annotation <- left_join(annotation, adduct_weights, by = "adduct")

  isotopes <- semi_join(annotation, adduct_weights, by = "adduct")
  isotopes <- purrr::pmap_dfr(isotopes, ~tibble())
  isotopes <- distinct(isotopes)

  annotation <- bind_rows(annotation, isotopes)
  annotation <- group_by(annotation, .data$compound)
  annotation <- mutate(annotation, score = 10^max(-Inf, weight, na.rm = TRUE))
  annotation <- ungroup(annotation)
  annotation <- compute_multiple_matches(annotation)
  annotation
}
