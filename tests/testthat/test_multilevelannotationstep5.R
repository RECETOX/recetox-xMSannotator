patrick::with_parameters_test_that(
  "multilevelannotationstep5",
  {
    if (exists("skip_function") && is.function(skip_function)) {
      skip_function()
    }
    
    testdata_dir <- file.path("test-data", subfolder)
    load(file.path(testdata_dir, "tempobjects.Rda"))

    testthat_wd <- getwd()

    outloc <- file.path(tempdir(), "multilevelannotationstep5", subfolder)
    dir.create(outloc, recursive = TRUE)
    file.copy(
      file.path(testdata_dir, "Stage4.csv"),
      file.path(outloc, "Stage4.csv")
    )
    expected <- read.csv(file.path(testdata_dir, "Stage5.csv"))


    multilevelannotationstep5(
      outloc = outloc,
      adduct_weights = adduct_weights)

    actual <- read.csv(file.path(outloc, "Stage5.csv"))
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