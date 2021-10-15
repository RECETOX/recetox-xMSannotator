#' @import tidyr
#' @import dplyr
#' @import plyr
#' @importFrom magrittr %>%
#' @export
get_chemscore <- function(chemical_id,
                          annotation,
                          correlation_threshold,
                          peak_correlation_matrix,
                          time_tolerance = 10,
                          adduct_weights,
                          filter_by = c("M+H"),
                          MplusH_abundance_ratio_check = TRUE,
                          outlocorig) {
  setwd(outlocorig)

  outloc1 <- paste(outlocorig, "/stage2/", sep = "")
  suppressWarnings(dir.create(outloc1))
  setwd(outloc1)

  if (length(annotation$mz) < 1) stop("No mz data found!")

  result <- compute_chemical_score(
    annotation,
    adduct_weights,
    global_cor,
    corthresh,
    filter_by,
    time_tolerance,
    chemical_id,
    MplusH_abundance_ratio_check
  )

  setwd("..")

  return(result$chemical_score)
}