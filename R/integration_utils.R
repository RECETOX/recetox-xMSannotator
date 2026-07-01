#' Recombine adduct formula and mass number difference.
#' @param formula Adduct formula, e.g 'M+H'
#' @param mass_number_difference Difference from monoisotopic mass.
#' @return Combined formula, e.g 'M+H_[+1]'.
construct_adduct_formula <- function(formula, mass_number_difference) {
  suffix <- sprintf("[%+1.0f]", mass_number_difference)
  formula <- paste(formula, suffix, sep = "_")
  return(formula)
}

#' Reformat annotation table with isotopes from `main` to be compatible with `master's` chemical score computing.
#'@param annotation Annotation table after [compute_isotopes()] with columns 
#' ['mz', 'rt', 'rt_cluster', 'module', 'multiple_match', 'expected_mass',
#' 'molecular_formula', 'adduct', 'mass_number_difference', 'monoisotopic_mass',
#' 'mean_intensity', 'mass_defect']
#' @return Reformatted annotation table for use with the `master` branch.
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
reformat_annotation_table <- function(annotation) {
  has_name <- "name" %in% colnames(annotation)

  master_annotation <- tibble(.rows = nrow(annotation))

  master_annotation <- master_annotation %>%
    mutate(
      Module_RTclust = paste(annotation$module, annotation$rt_cluster, sep = "_"),
      mz = annotation$mz,
      time = annotation$rt,
      MatchCategory = case_when(
        is.na(annotation$multiple_match) ~ "-",
        annotation$multiple_match ~ "Multiple",
        !annotation$multiple_match ~ "Unique"
      ),
      theoretical.mz = annotation$expected_mass,
      chemical_ID = paste("Formula", annotation$compound, sep = "_"),
      Formula = annotation$molecular_formula,
      MonoisotopicMass = annotation$monoisotopic_mass,
      Adduct = if_else(annotation$mass_number_difference == 0, annotation$adduct, construct_adduct_formula(annotation$adduct, annotation$mass_number_difference)),
      ISgroup = "-", # dummy column due to high complexity of reproducing the values
      time.y = annotation$rt,
      mean_int_vec = annotation$mean_intensity,
      MD = as.numeric(sprintf("0.%1.0f", annotation$mass_defect))
    )

  # Conditionally add Name column if it exists in input
  if (has_name) {
    master_annotation$Name <- annotation$name
  }

  return(master_annotation)
}

#' Reformat the peak correlation table row and column names to fit with `master` branch.
#' @param peak_table Table from which to take the mz & rt values.
#' @param peak_correlation_matrix Peak correlation table.
#' @return Peak correlation table where `peak` indices are replaced by `mz_rt`.
#' @importFrom magrittr set_colnames set_rownames
#' @export
reformat_correlation_matrix <- function(peak_table, peak_correlation_matrix, truncate = FALSE) {
  mz <- if(truncate) round(peak_table[, "mz"], 5) else peak_table[, "mz"]
  rt <- if(truncate) round(peak_table[, "rt"], 1) else peak_table[, "rt"]
  mz_rt <- paste0(mz, "_", rt)
  global_cor <- round(peak_correlation_matrix[,], 2)

  global_cor <- magrittr::set_colnames(global_cor, mz_rt)
  global_cor <- magrittr::set_rownames(global_cor, mz_rt)
  return(global_cor)
}
