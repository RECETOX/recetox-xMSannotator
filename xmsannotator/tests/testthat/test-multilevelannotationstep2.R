patrick::with_parameters_test_that(
  "multilevelannotationstep2:",
  {
    # Arrange
    testthat_wd <- getwd()
    test_path <- file.path(
      testthat_wd,
      "test-data",
      test_identifier
    )

    expected <- readRDS(file.path(test_path, "chemscoremat.Rds"))
    load(file = file.path(test_path, "step1_results.Rda"))
    load(file = file.path(test_path, "global_cor.Rda"))
    load(file = file.path(test_path, "tempobjects.Rda"))

    outloc <- file.path(
      tempdir(),
      "multilevelannotationstep2",
      test_identifier
    )
    dir.create(outloc, recursive = TRUE)

    # Act
    actual <- lapply(
      1:num_sets,
      call_multilevelannotationstep2,
      outloc = outloc,
      max.rt.diff = max.rt.diff,
      chemids_split = chemids_split,
      num_sets = num_sets,
      mchemdata = mchemdata,
      mass_defect_window = mass_defect_window,
      corthresh = corthresh,
      global_cor = global_cor,
      mzid = mzid,
      adduct_table = adduct_table,
      adduct_weights = adduct_weights,
      max_isp = max_isp,
      MplusH.abundance.ratio.check = MplusH.abundance.ratio.check,
      mass_defect_mode = mass_defect_mode,
      chemids = chemids,
      isop_res_md = isop_res_md,
      filter.by = filter.by
    )
    actual <- plyr::ldply(actual, rbind)
    actual$Formula <- gsub(actual$Formula, pattern = "_.*", replacement = "")


    keys <- c("mz", "time", "Name", "Adduct", "Formula", "chemical_ID", "cur_chem_score")

    actual$MonoisotopicMass <- as.character(actual$MonoisotopicMass)
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
      mismatchCount = 1000
    )

    # Annihilate
    setwd(outloc)
    unlink("stage2", recursive = TRUE)
    setwd(testthat_wd)

    # Assert
    expect_equal(actual, expected)

    rm(actual)
    gc(reset = TRUE)
  },
  patrick::cases(
    qc_solvent = list(test_identifier = "qc_solvent"),
    batch1_neg = list(test_identifier = "batch1_neg"),
    sourceforge = list(test_identifier = "sourceforge"),
    qc_matrix = list(test_identifier = "qc_matrix")
  )
)