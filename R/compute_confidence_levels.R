#' @import dplyr
#' @importFrom rlang .data
compute_expected_confidences <- function (annotation, expected_adducts) {
  annotation <- mutate(annotation, expected_adduct = .data$adduct %in% {{ expected_adducts }})
  annotation <- group_by(annotation, .data$compound)
  annotation <- mutate(annotation, confidence = ifelse(any(.data$expected_adduct), 2, 0))
  annotation <- ungroup(annotation)
  annotation
}

#' Boost confidence level of certain annotations
#'
#' For certain annotations set the confidence level to the highest level (4). The annotations to be boosted are selected
#' based on the compound ID match. The selection can be further narrowed by specifing the optional \emph{mz} and
#' \emph{rt} values in the \emph{boosted_compounds} parameter.
#'
#' @param annotation data frame of annotations required to have \emph{compound}, \emph{mz}, \emph{rt},
#' \emph{confidence}, and \emph{score} columns
#' @param boosted_compounds data frame of compounds to be boosted required to have 3 columns \emph{compound}, \emph{mz},
#' \emph{rt}. The \emph{mz} and \emph{rt} columns are allowed to contain NAs causing the boosting to be less restrictive
#' @param mass_tolerance how close the \emph{mz} values must be to each other in order to pronounce them as a match
#' @param time_tolerance how close the \emph{rt} values must be to each other in order to pronounce them as a match
#'
#' @return annotation with increased confidence levels and scores for boosted compounds and with a new added column
#' \emph{boosted} indicating whether the annotation entry was boosted or not
#'
#' @import dplyr
#' @import tidyr
#' @importFrom rlang .data
compute_boosted_confidences <- function(annotation, boosted_compounds, mass_tolerance, time_tolerance) {
  is_rt_near <- function (x, y) near(x, y, tol = time_tolerance)
  is_mz_near <- function (x, y) near(x, y, tol = max(abs(x), abs(y)) * mass_tolerance)

  boosted <- sapply(seq_len(nrow(annotation)), function (i) {
    id_mask <- annotation$compound[[i]] == boosted_compounds$compound
    rt_mask <- replace_na(is_rt_near(annotation$rt[[i]], boosted_compounds$rt), TRUE)
    mz_mask <- replace_na(is_mz_near(annotation$mz[[i]], boosted_compounds$mz), TRUE)
    any(id_mask & rt_mask & mz_mask)
  })

  annotation <- mutate(annotation, boosted = {{ boosted }})
  annotation <- mutate(annotation, confidence = ifelse(.data$boosted, 4, .data$confidence))
  annotation <- mutate(annotation, score = ifelse(.data$boosted, 100 * .data$score, .data$score))
  annotation
}


compute_confidence_levels <- function(annotation, expected_adducts, boosted_compounds, mass_tolerance, time_tolerance) {
  # monoisotopic_adduct <- tibble::tibble(adduct = "M", charge = 1, mass = 0, molecules = 1)
  # adduct_table <- dplyr::bind_rows(adduct_table, monoisotopic_adduct)

  # dplyr::group_by(compound) %>%
  # dplyr::mutate(confidence = if (n() > 1 & all(score < 10) & n_distinct(adducts_with_isotopes) < 2) 0 else confidence) %>%
  # dplyr::ungroup()

  if (!is.null(expected_adducts)) {
    annotation <- compute_expected_confidences(annotation, expected_adducts)
  }

  if (!is.null(boosted_compounds)) {
    annotation <- compute_boosted_confidences(annotation, boosted_compounds, mass_tolerance, time_tolerance)
  }

  annotation <- compute_multiple_matches(annotation)
  annotation
}
