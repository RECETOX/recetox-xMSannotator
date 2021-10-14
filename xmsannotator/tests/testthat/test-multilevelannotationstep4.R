patrick::with_parameters_test_that(
  "multilevelannotation step 4 works",
  {
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
    qc_matrix = list(subfolder = "qc_matrix"),
    batch1_neg = list(subfolder = "batch1_neg"),
    sourceforge = list(subfolder = "sourceforge")
  )
)