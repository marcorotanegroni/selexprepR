test_that("round inference follows the conservative metadata cascade", {
    runs <- S4Vectors::DataFrame(
        run_accession = c("SRR_attr", "SRR_title", "SRR_library", "SRR_glued", "SRR_ambiguous"),
        sample_title = c("opaque", "Thrombin Round 2", "opaque", "DNAFOXR00", "Round 1 Round 2"),
        library_name = c("", "", "RAPT26-7R", "", ""),
        sample_attributes = I(list(
            c(selection_round = "4"), NULL, NULL, NULL, NULL
        ))
    )

    assignments <- infer_selexprep_rounds(runs)

    expect_equal(assignments$round_number, c(4L, 2L, 7L, 0L, 1L))
    expect_identical(assignments$confidence, c("HIGH", "HIGH", "MEDIUM", "HIGH", "MEDIUM"))
    expect_identical(assignments$source_field[[1L]], "sample_attributes")
    expect_identical(assignments$target_hint[[2L]], "Thrombin")
    expect_identical(assignments$round_candidates[[5L]], c(1L, 2L))
    expect_false(selexprepR:::.safe_round_assignments(assignments)[[5L]])
})

test_that("round overrides are explicit and unknown metadata remains unresolved", {
    runs <- data.frame(
        run_accession = c("SRR_manual", "SRR_unknown"),
        sample_title = c("no useful annotation", "still opaque")
    )

    assignments <- infer_selexprep_rounds(runs, overrides = c(SRR_manual = 9L))

    expect_identical(assignments$round_number, c(9L, NA_integer_))
    expect_identical(assignments$source_field[[1L]], "seed_override")
    expect_identical(assignments$confidence[[2L]], "NONE")
    expect_error(infer_selexprep_rounds(runs, overrides = c(3L)), "named vector")
})
