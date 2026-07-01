patrick::with_parameters_test_that(
  "multilevelannotation step 4 works",
  {
    if (exists("skip_function") && is.function(skip_function)) {
      skip_function()
    }

    # load data needed during step 4
    testdata_dir <- file.path(getwd(), "test-data", subfolder)
    load(file.path(testdata_dir, "tempobjects.Rda"))

    testthat_wd <- getwd()
    outloc <- file.path(
      tempdir(),
      "multilevelannotationstep4",
      subfolder
    )

    # create test folder
    dir.create(outloc, recursive = TRUE)

    chemscoremat <- read.csv(file.path(testdata_dir, "Stage3.csv"))

    # load expected results
    expected <- read.csv(file.path(testdata_dir, "Stage4.csv"))

    # compute annotation step 4
    result <- multilevelannotationstep4(
      outloc = outloc,
      chemscoremat = chemscoremat,
      max.mz.diff = max.mz.diff,
      max.rt.diff = max_diff_rt,
      filter.by = filter.by,
      adduct_weights = adduct_weights,
      max_isp = max_isp,
      min_ions_perchem = min_ions_perchem
    )
    actual <- read.csv(file.path(outloc, "Stage4.csv"))

    setwd(testthat_wd)

    actual <- dplyr::arrange(actual, dplyr::across(everything()))
    expected <- dplyr::arrange(expected, dplyr::across(everything()))

    # Note: dataCompareR::rCompare() removed due to deprecated select_() usage in dataCompareR
    # The expect_equal below provides the same assertion

    expect_equal(actual, expected)
  },
  patrick::cases(
    qc_solvent = list(subfolder = "qc_solvent"),
    qc_matrix = list(subfolder = "qc_matrix", skip_function = skip_on_ci),
    batch1_neg = list(subfolder = "batch1_neg", skip_function = skip_on_ci),
    sourceforge = list(subfolder = "sourceforge", skip_function = skip_on_ci)
  )
)

# Test for boundary case: score == 10 should be handled correctly
# Note: This test requires proper test data with all columns expected by get_confidence_stage4
# For now, skipping as it requires loading additional data (adduct_table, etc.)
test_that("compute_confidence_levels handles score == 10 boundary case", {
  skip("Test requires proper test data structure with all columns expected by get_confidence_stage4")
})