patrick::with_parameters_test_that("basic advanced_annotation functionality", {
  if (exists("skip_function") && is.function(skip_function)) {
    skip_function()
  }
  skip("")

  testthat_wd <- getwd()
  outloc <- file.path(tempdir(), "advanced_annotation")
  dir.create(outloc, recursive = TRUE)

  expected <- read.csv(expected_output)

  peaks <- arrow::read_parquet(peak_table)
  DB <- arrow::read_parquet(database)
  adduct_table <- sample_adduct_table %>% filter(adduct %in% adduct_names)

  annotation <- advanced_annotation(
    peaks,
    DB,
    adduct_table = adduct_table,
    peak_rt_width = peak_rt_width,
    time_tolerance = time_tolerance,
    intensity_deviation_tolerance = intensity_deviation_tolerance,
    outloc = outloc
  )

  actual <- read.csv(file.path(outloc, "Stage5.csv"))
  setwd(testthat_wd)

  expect_equal(actual, expected)
},
  patrick::cases(
    simple_case = list(test_identifier = "simple_case",
                       expected_output = "test-data/advanced_annotation/expected_small.csv",
                       peak_table = "test-data/advanced_annotation/peak_table.parquet",
                       database = "test-data/advanced_annotation/database_small.parquet",
                       adduct_names = c("M-H", "M+H"),
                       peak_rt_width = 1,
                       time_tolerance = 10,
                       intensity_deviation_tolerance = 0.1),
    large_case = list(test_identifier = "large_case",
                      expected_output = "test-data/advanced_annotation/expected.csv",
                      peak_table = "test-data/batch1_neg.parquet",
                      database = "test-data/advanced_annotation/hmdb.parquet",
                      adduct_names = c("M-H", "2M-H", "M-2H"),
                      peak_rt_width = 15,
                      time_tolerance = 15,
                      intensity_deviation_tolerance = 0.2,
                      skip_function = skip_on_ci)
  )
)
