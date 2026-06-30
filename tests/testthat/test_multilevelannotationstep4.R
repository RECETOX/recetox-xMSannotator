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

    actual <- dplyr::arrange_all(actual)
    expected <- dplyr::arrange_all(expected)

    comparison <- dataCompareR::rCompare(
      actual,
      expected,
      keys = names(actual)
    )

    dataCompareR::saveReport(
      comparison,
      reportName = subfolder,
      reportLocation = outloc,
      showInViewer = FALSE,
      mismatchCount = 1000
    )

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
test_that("compute_confidence_levels handles score == 10 boundary case", {
  # Create minimal test data with score exactly equal to 10
  curdata_score_10 <- tibble::tibble(
    chemical_ID = "TEST001",
    Adduct = "[M+H]+",
    mz = "123.456",
    time = 100,
    theoretical.mz = "123.450",
    score = 10,  # Exactly 10 - boundary case
    Formula = "C6H12O6"
  )

  # Create test data with score just below 10
  curdata_score_9 <- tibble::tibble(
    chemical_ID = "TEST002",
    Adduct = "[M+H]+",
    mz = "123.456",
    time = 100,
    theoretical.mz = "123.450",
    score = 9,  # Just below 10
    Formula = "C6H12O6"
  )

  # Create test data with score just above 10
  curdata_score_11 <- tibble::tibble(
    chemical_ID = "TEST003",
    Adduct = "[M+H]+",
    mz = "123.456",
    time = 100,
    theoretical.mz = "123.450",
    score = 11,  # Just above 10
    Formula = "C6H12O6"
  )

  # Create minimal adduct_weights
  adduct_weights <- tibble::tibble(
    V1 = c("[M+H]+", "[M+Na]+"),
    V2 = c(1, 0)
  )

  # Test that score == 10 is processed without error
  result_10 <- compute_confidence_levels(
    c = 1,
    chemids = "TEST001",
    chemscoremat = as.data.frame(curdata_score_10),
    filter.by = NA,
    max.rt.diff = 30,
    adduct_weights = adduct_weights,
    max_isp = 5,
    min_ions_perchem = 1
  )

  # Test that score < 10 is processed without error
  result_9 <- compute_confidence_levels(
    c = 1,
    chemids = "TEST002",
    chemscoremat = as.data.frame(curdata_score_9),
    filter.by = NA,
    max.rt.diff = 30,
    adduct_weights = adduct_weights,
    max_isp = 5,
    min_ions_perchem = 1
  )

  # Test that score > 10 is processed without error
  result_11 <- compute_confidence_levels(
    c = 1,
    chemids = "TEST003",
    chemscoremat = as.data.frame(curdata_score_11),
    filter.by = NA,
    max.rt.diff = 30,
    adduct_weights = adduct_weights,
    max_isp = 5,
    min_ions_perchem = 1
  )

  # All results should have Confidence column and not error
  expect_true("Confidence" %in% colnames(result_10))
  expect_true("Confidence" %in% colnames(result_9))
  expect_true("Confidence" %in% colnames(result_11))

  # All results should have chemical_ID column
  expect_true("chemical_ID" %in% colnames(result_10))
  expect_true("chemical_ID" %in% colnames(result_9))
  expect_true("chemical_ID" %in% colnames(result_11))
})