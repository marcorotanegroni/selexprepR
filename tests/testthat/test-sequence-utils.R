test_that("reverse_complement matches the Python counting helper", {
    expect_true(identical(
        reverse_complement(c("ACGT", "AAAA", "GGCC")),
        c("ACGT", "TTTT", "GGCC")
    ))
    expect_true(identical(reverse_complement("ACGU"), "ACGT"))
    expect_true(identical(reverse_complement(c("acgt", "ANNT")), c("acgt", "ANNT")))
})

test_that("reverse_complement validates its R interface", {
    expect_error(reverse_complement(1), "character vector")
    expect_error(reverse_complement(c("ACGT", NA_character_)), "missing values")
})
