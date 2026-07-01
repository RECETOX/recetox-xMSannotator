test_that("simple_annotation functionality on sample data", {
  peaks <- arrow::read_parquet("test-data/simple_annotation_sample/peaks.parquet")
  DB <- arrow::read_parquet("test-data/simple_annotation_sample/compounds.parquet")
  expected <- arrow::read_parquet("test-data/simple_annotation_sample/expected.parquet")

  result <- simple_annotation(peaks, DB)

  expect_equal(result, expected)
})

# Test that check_element correctly handles element symbols that are prefixes of others
# e.g., "C" should not match in "Cl", "Ca", "Co", etc.
test_that("check_element handles element prefix edge cases correctly", {
  # Carbon should not be detected in chlorine, calcium, or cobalt
  expect_equal(check_element("Cl", "C"), 0)
  expect_equal(check_element("Cl2", "C"), 0)
  expect_equal(check_element("Ca", "C"), 0)
  expect_equal(check_element("Co", "C"), 0)
  expect_equal(check_element("NaCl", "C"), 0)

  # But carbon should be correctly counted when present
  expect_equal(check_element("C6H12O6", "C"), 6)
  expect_equal(check_element("CCl4", "C"), 1)
  expect_equal(check_element("C2H5Cl", "C"), 2)
  expect_equal(check_element("CH3COOH", "C"), 2)
  expect_equal(check_element("Na2CO3", "C"), 1)
  expect_equal(check_element("CaCO3", "C"), 1)

  # Similar tests for other elements that could have prefix issues
  expect_equal(check_element("Na", "N"), 0)  # N should not match in Na
  expect_equal(check_element("Si", "S"), 0)  # S should not match in Si
  expect_equal(check_element("Pb", "P"), 0)  # P should not match in Pb

  # Correct counts for these elements
  expect_equal(check_element("NaCl", "Na"), 1)
  expect_equal(check_element("H2SiO3", "Si"), 1)
  expect_equal(check_element("Pb(NO3)2", "Pb"), 1)
})
