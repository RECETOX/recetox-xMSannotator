test_that("Integration utils: annotation table reformating", {
  main_annotation_table <- readRDS("test-data/integration_utils/main_annotation_table.rds")

  actual <- reformat_annotation_table(main_annotation_table)
  expected <- readRDS("test-data/integration_utils/master_annotation_table.rds")

  # Use across() instead of deprecated arrange_all()
  actual <- dplyr::arrange(actual, dplyr::across(everything()))
  expected <- dplyr::arrange(expected, dplyr::across(everything()))

  # Note: dataCompareR::rCompare() removed due to deprecated select_() usage in dataCompareR
  # The expect_equal below provides the same assertion
  expect_equal(actual, expected)
})

patrick::with_parameters_test_that("Global cor adapter works", {
  if (exists("skip_function") && is.function(skip_function)) {
    skip_function()
  }
  
  testthat_wd <- getwd()
  test_path <- file.path(testthat_wd, "test-data")

  load(file.path(test_path, test_identifier, "global_cor_integ.Rda"))

  peaks_filename <- file.path(test_path, paste0(test_identifier, ".parquet"))
  peak_table <- load_peak_table_parquet(peaks_filename)
  peak_intensity_matrix <- get_peak_intensity_matrix(peak_table)
  peak_correlation_matrix <- compute_peak_correlations(peak_intensity_matrix, correlation_method = "p")

  actual <- reformat_correlation_matrix(peak_table, peak_correlation_matrix, truncate = TRUE)

  expect_true(all.equal(actual, global_cor, check.names = TRUE, countEQ = TRUE))
},
  patrick::cases(
    qc_solvent = list(test_identifier = "qc_solvent"),
    batch1_neg = list(test_identifier = "batch1_neg", skip_function = skip_on_ci),
    sourceforge = list(test_identifier = "sourceforge", skip_function = skip_on_ci),
    qc_matrix = list(test_identifier = "qc_matrix", skip_function = skip_on_ci)
  )
)
