test_that("simple_annotation", {
  df <- simple_annotation(sample_peak_table, sample_compound_table)

  expect_s3_class(df, 'data.frame')
})
