patrick::with_parameters_test_that("Compute chemscore can be called isolated", {
    # Arrange
    testthat_wd <- getwd()
    outloc <- file.path(tempdir(), "get_chemscore", test_identifier)
    test_path <- file.path(
      testthat_wd,
      "test-data",
      test_identifier
    )

    expected <- readRDS(file.path(test_path, "chemscoremat.Rds"))
    if(!dir.exists(outloc)) {
        dir.create(outloc, recursive = TRUE)
    }

    annotation_file <- file.path("test-data", "get_chemscore",paste0(test_identifier, "_annotation.Rds"))
    isotopes <- readRDS(annotation_file)  

    if(is.factor(isotopes$MonoisotopicMass)) {
      isotopes$MonoisotopicMass <- as.numeric(levels(isotopes$MonoisotopicMass))[isotopes$MonoisotopicMass]
    }

    isotopes$MonoisotopicMass <- as.numeric(isotopes$MonoisotopicMass)
    isotopes$theoretical.mz <- as.numeric(isotopes$theoretical.mz)


    data(adduct_weights)
    load(file = file.path(test_path, "global_cor.Rda"))

    annotation <- isotopes  

    actual <- purrr::pmap_dfr(
        annotation,
        ~ get_chemscore(...,
            annotation = isotopes,
            adduct_weights = adduct_weights,
            correlation_threshold = 0.7,
            peak_correlation_matrix = global_cor,
            time_tolerance = max_diff_rt,
            outlocorig = outloc
      )
    )

    actual$Formula <- gsub(actual$Formula, pattern = "_.*", replacement = "")
    keys <- colnames(actual)

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

    comparison <- dataCompareR::rCompare(actual, expected, keys = keys)
    dataCompareR::saveReport(
      comparison,
      reportName = test_identifier,
      reportLocation = outloc,
      showInViewer = FALSE,
      mismatchCount = 10000
    )

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
    batch1_neg = list(test_identifier = "batch1_neg", max_diff_rt = 0.5),
    sourceforge = list(test_identifier = "sourceforge", max_diff_rt = 2),
    qc_matrix = list(test_identifier = "qc_matrix", max_diff_rt = 0.5)
  )
)
