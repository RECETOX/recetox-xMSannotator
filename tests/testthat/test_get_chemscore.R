patrick::with_parameters_test_that("Compute chemscore can be called isolated", {

    if (exists("skip_function") && is.function(skip_function)) {
      skip_function()
    }

    # Arrange
    testthat_wd <- getwd()
    test_path <- file.path(
      testthat_wd,
      "test-data",
      test_identifier
    )
    browser()
    load(file = file.path(test_path, "tempobjects.Rda"))

    outloc <- file.path(tempdir(), "get_chemscore", test_identifier)

    if(!dir.exists(outloc)) {
      dir.create(outloc, recursive = TRUE)
    }

    expected <- readRDS(file.path(test_path, "expected_chemscore.Rds"))

    setwd(testthat_wd)

    annotation_file <- file.path("test-data", "get_chemscore", paste0(test_identifier, "_annotation.Rds"))
    isotopes <- readRDS(annotation_file)

    if(is.factor(isotopes$MonoisotopicMass)) {
      isotopes$MonoisotopicMass <- as.numeric(levels(isotopes$MonoisotopicMass))[isotopes$MonoisotopicMass]
    }

    load(file = file.path(test_path, "global_cor.Rda"))

    actual <- purrr::pmap_dfr(
        isotopes,
        ~get_chemscore(...,
          annotation = isotopes,
          mass_defect_window = mass_defect_window,
          adduct_weights = adduct_weights,
          corthresh = corthresh,
          global_cor = global_cor,
          max_diff_rt = max_diff_rt,
          outlocorig = outloc,
          filter.by = c("M-H", "M+H")
        )
    )
    actual <- unique(actual)

    keys <- c("mz", "time", "Name", "Adduct", "Formula", "chemical_ID", "cur_chem_score", "MatchCategory")

    actual$MonoisotopicMass[is.na(actual$MonoisotopicMass)] <- "-"
    actual$theoretical.mz[is.na(actual$theoretical.mz)] <- "-"

    actual$MonoisotopicMass <- as.character(actual$MonoisotopicMass)
    actual$theoretical.mz <- as.character(actual$theoretical.mz)
    expected$Module_RTclust <- as.character(expected$Module_RTclust)

    actual <- dplyr::arrange_at(
      actual, keys
    )
    expected <- dplyr::arrange_at(
      expected, keys
    )

    # comparison <- dataCompareR::rCompare(actual, expected, keys = keys)
    # dataCompareR::saveReport(
    #   comparison,
    #   reportName = test_identifier,
    #   reportLocation = outloc,
    #   showInViewer = FALSE,
    #   mismatchCount = 10000
    # )

    write.csv(actual, file = file.path(outloc, "chemscoremat_actual.csv"))
    write.csv(expected, file = file.path(outloc, "chemscoremat_expected.csv"))

    setwd(testthat_wd)

    # Assert
    expect_equal(actual, expected)

    rm(actual)
    gc(reset = TRUE)
},
  patrick::cases(
    qc_solvent = list(test_identifier = "qc_solvent", max_diff_rt = 0.5),
    batch1_neg = list(test_identifier = "batch1_neg", max_diff_rt = 0.5, skip_function = skip),
    sourceforge = list(test_identifier = "sourceforge", max_diff_rt = 2, skip_function = skip_on_ci),
    qc_matrix = list(test_identifier = "qc_matrix", max_diff_rt = 0.5, skip_function = skip_on_ci)
  )
)


test_that("Chemscores range from 0 to Inf", {
  skip_on_ci()
  test_identifier <- "batch1_neg"

  testthat_wd <- getwd()
  test_path <- file.path(
    testthat_wd,
    "test-data",
    test_identifier
  )

  outloc <- file.path(tempdir(), "get_chemscore_2", test_identifier)

  if(!dir.exists(outloc)) {
    dir.create(outloc, recursive = TRUE)
  }

  annotation <- readRDS(file.path("test-data", "get_chemscore", "annotation_with_isotopes_v3.Rds"))
  global_cor <- as.data.frame(readRDS(file.path("test-data", "get_chemscore", "global_cor.Rds")))

  adduct_names <- c("M-H", "2M-H")
  adduct_weights <- as.data.frame(tibble(adduct = adduct_names, weight = rep_len(5, length(adduct_names))))

  actual <- purrr::pmap_dfr(
    annotation,
    ~ get_chemscore(...,
      annotation = annotation,
      mass_defect_window = 0.01,
      adduct_weights = adduct_weights,
      corthresh = 0.7,
      global_cor = global_cor,
      max_diff_rt = 15,
      outlocorig = outloc,
      filter.by = c("M-H")
    )
  )
  expect_equal(min(actual$cur_chem_score), 0)
  expect_gt(max(actual$cur_chem_score), 10^5)
})