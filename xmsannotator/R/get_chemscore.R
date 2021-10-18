#' @import tidyr
#' @import dplyr
#' @import plyr
#' @importFrom magrittr %>%
#' @export
get_chemscore <- function(...,
                          annotation,
                          peak_correlation_matrix,
                          adduct_weights,
                          correlation_threshold,
                          time_tolerance = 10,
                          filter_by = c("M+H"),
                          MplusH_abundance_ratio_check = TRUE,
                          outlocorig) {
  setwd(outlocorig)

  query <- tibble(...)
  outloc1 <- paste(outlocorig, "/stage2/", sep = "")
  suppressWarnings(dir.create(outloc1))
  setwd(outloc1)

  if (length(annotation$mz) < 1) stop("No mz data found!")

  mchemicaldata <- annotation %>% filter(
    .$chemical_ID == query$chemical_ID
  )

  if(nrow(mchemicaldata %>% filter(is.na(MonoisotopicMass))) > 0) {
    mchemicaldata$MonoisotopicMass[is.na(mchemicaldata$MonoisotopicMass)] <- "-"
    mchemicaldata$theoretical.mz[is.na(mchemicaldata$theoretical.mz)] <- "-"

    mchemicaldata$MonoisotopicMass <- as.character(mchemicaldata$MonoisotopicMass)
    mchemicaldata$theoretical.mz <- as.character(mchemicaldata$theoretical.mz)
  }

  result <- compute_chemical_score(
    mchemicaldata,
    adduct_weights,
    peak_correlation_matrix,
    correlation_threshold,
    filter_by,
    time_tolerance,
    query$chemical_ID,
    MplusH_abundance_ratio_check
  )

  setwd("..")

  query$cur_chem_score <- result$chemical_score
  return(query)
}