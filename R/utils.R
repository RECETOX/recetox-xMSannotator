lexicographic_rank <- function(...) {
  .o <- order(...)
  .x <- cbind(...)[.o,]
  cumsum(!duplicated(.x))[order(.o)]
}

#' @import dplyr
as_peak_table <- function(data, intensities = FALSE) {
  optional <- rlang::quo(any_of('peak'))
  required <- rlang::quo(any_of(c('mz', 'rt')))
  additional <- rlang::quo(starts_with('intensity'))
  data <- select(data, !!optional, !!required, if (intensities) !!additional else NULL)

  if (!is.element('peak', colnames(data))) {
    data$peak <- lexicographic_rank(data$mz, data$rt)
  }

  stopifnot(anyDuplicated(data$peak) == 0)
  stopifnot(is.integer(data$peak))
  stopifnot(is.numeric(data$mz))
  stopifnot(is.numeric(data$rt))

  return(data)
}

#' @import dplyr
as_adduct_table <- function(data) {
  data <- select(data, 'adduct', 'charge', 'factor', 'mass')

  stopifnot(anyDuplicated(data$adduct) == 0)
  stopifnot(is.integer(data$charge))
  stopifnot(is.integer(data$factor))
  stopifnot(is.numeric(data$mass))

  return(data)
}

#' @import dplyr
as_compound_table <- function(data) {
  optional <- c('compound', 'recetox_cid')
  required <- c('monoisotopic_mass', 'molecular_formula')
  data <- select(data, any_of(optional), all_of(required))

  if ('recetox_cid' %in% colnames(data)) {
    data <- mutate(data, compound = .data$recetox_cid)
  }

  stopifnot('compound' %in% colnames(data))
  stopifnot(anyDuplicated(data$compound) == 0)
  stopifnot(is.numeric(data$monoisotopic_mass))
  stopifnot(is.character(data$molecular_formula))

  return(data)
}

#' @import dplyr
as_expected_adducts_table <- function (data) {
  data <- select(data, 'adduct')
  return(data)
}

#' @import dplyr
as_boosted_compounds_table <- function (data) {
  data <- select(data, all_of('compound'), any_of(c('mz', 'rt')))

  stopifnot(is.numeric(data$mz))
  stopifnot(is.numeric(data$rt))

  return(data)
}

#' @import dplyr
load_parquet <- function (file, columns) {
  rlang::with_handlers(
    arrow::read_parquet(file, col_select = any_of(columns)),
    error = ~ rlang::abort(paste("The file", toString(file), "seams not to be a valid Parquet file."), parent = .)
  )
}

load_csv <- function (file, columns) {
  data <- readr::read_csv(file)
  data <- select(data, any_of(columns))
}

#' @export
load_peak_table_parquet <- function(file) {
  data <- load_parquet(file, columns = c("peak", "mz", "rt"))
  as_peak_table(data)
}

#' @export
load_adduct_table_parquet <- function(file) {
  data <- load_parquet(file, columns = c("adduct", "charge", "mass", "n_molecules"))
  as_adduct_table(data)
}

#' @export
load_compound_table_parquet <- function(file) {
  data <- load_parquet(file, columns = c("compound", "recetox_cid", "monoisotopic_mass", "molecular_formula"))
  as_compound_table(data)
}

#' @export
load_expected_adducts_csv <- function (file) {
  data <- load_csv(file, columns = "adduct")
  as_expected_adducts_table(data)
}

#' @export
load_boost_compounds_csv <- function (file) {
  data <- load_csv(file, columns = c("compound", "mz", "rt"))
  as_boosted_compounds_table(data)
}

#' @export
load_peak_table_hdf <- function(file, intensities = FALSE) {
  data <- rhdf5::h5read(file, "peaks")
  as_peak_table(data, intensities = intensities)
}

#' @export
save_parquet <- function(data, file) {
  invisible(arrow::write_parquet(data, file))
}
