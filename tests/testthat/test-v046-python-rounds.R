test_that("R round assignments match Python v0.4.6 on 41 identical metadata inputs", {
    reference <- jsonlite::fromJSON(test_path("fixtures", "v046-round-reference.json"),
        simplifyVector = FALSE)
    expect_identical(reference$source$commit, "63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5")
    expect_length(reference$cases, 41L)
    for (index in seq_along(reference$cases)) {
        case <- reference$cases[[index]]
        arguments <- c(list(run_accession = "SRR1"), case$input)
        observed <- do.call(selexprepR:::.single_round_assignment, arguments)
        expected <- case$expected
        number <- if (is.null(expected$round_number)) NA_integer_ else as.integer(expected$round_number)
        candidates <- as.integer(unlist(expected$round_candidates, use.names = FALSE))
        expect_identical(observed$round_number, number, info = as.character(index))
        expect_identical(observed$confidence, expected$confidence, info = as.character(index))
        expect_identical(observed$source_field, expected$source_field, info = as.character(index))
        expect_identical(observed$candidates, candidates, info = as.character(index))
        expect_identical(is.na(observed$round_number) || length(observed$candidates) != 1L,
            expected$is_unassigned, info = as.character(index))
    }
})
