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
    annotation <- readRDS(annotation_file)
    # annotation$peak <- as.integer(rownames(annotation))


    peaks_filename <- file.path("test-data/get_chemscore", paste0(test_identifier,".rds"))

    data(adduct_weights)
    peaks <- readRDS(peaks_filename)
    peaks$peak <- rownames(peaks)
    peaks <- dplyr::rename(peaks, rt = time)
    peaks <- dplyr::rename_with(
        peaks,
        ~ paste0("intensity_", .x),
        starts_with(c("Tribrid", "Thermo", "rep"))
    )
    peaks$peak <- as.integer(peaks$peak)

    intensity_matrix <- t(select(peaks, starts_with("intensity")))
    intensity_matrix <- magrittr::set_colnames(intensity_matrix, paste0(round(peaks$mz, 5), "_", round(peaks$rt, 1)))
    correlation_matrix <- WGCNA::cor(intensity_matrix, use = "p", method = "p")

    chemical_groups <- group_split(annotation, chemical_ID)
    chemical_groups_v2 <- group_by(annotation, chemical_ID)

    # actual <- annotation %>% group_by(chemical_ID) %>%
    #     mutate(
    #     cur_chem_score = get_chemscore(
    #         chemical_ID,
    #         annotation = annotation,
    #         adduct_weights = adduct_weights,
    #         correlation_threshold = 0.7,
    #         peak_correlation_matrix = correlation_matrix,
    #         time_tolerance = 10,
    #         outlocorig = outloc
    #     )
    # )
    chemical_IDs <- as.list(annotation$chemical_ID)

    actual <- purrr::pmap_dfr(
        annotation,
        ~ get_chemscore(...,
            annotation = annotation,
            adduct_weights = adduct_weights,
            correlation_threshold = 0.7,
            peak_correlation_matrix = correlation_matrix,
            time_tolerance = 10,
            outlocorig = outloc
        )
    )

    write.csv(actual, file = file.path(outloc, "chemscoremat.csv"))

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

    setwd(testthat_wd)

    # Assert
    expect_equal(actual, expected)

    rm(actual)
    gc(reset = TRUE)
},
  patrick::cases(
    qc_solvent = list(test_identifier = "qc_solvent"),
    batch1_neg = list(test_identifier = "batch1_neg"),
    #sourceforge = list(test_identifier = "sourceforge"),
    qc_matrix = list(test_identifier = "qc_matrix")
  )
)
