test_that("R extracts the same sequences as Python v0.4.6 on identical tagged reads", {
    references <- jsonlite::fromJSON(
        test_path("fixtures", "v046-extraction-reference.json"),
        simplifyVector = FALSE
    )
    for (name in names(references)) {
        reference <- references[[name]]
        fields <- reference$report
        as_named_numeric <- function(values) {
            stats::setNames(as.numeric(unlist(values)), names(values))
        }
        fields["known_adapter_hits"] <- list(as_named_numeric(fields$known_adapter_hits))
        fields["n_length_distribution"] <- list(as_named_numeric(fields$n_length_distribution))
        report <- selexprepR:::.new_library_report(fields)
        pools <- lapply(reference$sequences_by_round, unlist, use.names = FALSE)
        extraction <- selexprep_extract(pools, report)
        observed <- lapply(extraction$sequences_by_round, as.character)
        expected <- lapply(reference$expected_sequences_by_round, unlist, use.names = FALSE)
        expect_identical(observed, expected, info = name)
        expect_length(extraction$low_yield_rounds, 0L)
    }
})
