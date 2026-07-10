.report_with_mode <- function(mode, primer_5p = "GGTAATACGACTCACTATAGGG",
    primer_3p = "CCATGCATGCATGCATGCAT", orientation = "FORWARD") {
    fields <- unclass(read_library_report(test_path("fixtures", "library-report-both-primers.json")))
    fields["primer_5p"] <- list(primer_5p)
    fields["primer_3p"] <- list(primer_3p)
    fields$orientation <- orientation
    fields$extraction_mode <- mode
    fields$full_insert_recovered <- identical(mode, "BOTH_PRIMERS_SINGLE_READ")
    fields$required_action <- switch(mode,
        BOTH_PRIMERS_SINGLE_READ = "NONE",
        FIVE_PRIME_ONLY = "NONE",
        THREE_PRIME_ONLY = "NONE",
        PAIRED_END_SPLIT_PRIMERS = "READ_MERGING_RECOMMENDED"
    )
    selexprep:::.new_library_report(fields)
}

test_that("selexprep_extract retains only fully linked-primer reads", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    good <- paste0(primer_5p, "ACGTACGT", primer_3p)
    one_mismatch <- paste0("A", substr(primer_5p, 2L, nchar(primer_5p)), "TTTT", primer_3p)
    bad <- paste0(primer_5p, "GGGGGGGG")

    extraction <- selexprep_extract(
        list(round_00 = c(good, one_mismatch, bad)),
        .report_with_mode("BOTH_PRIMERS_SINGLE_READ")
    )

    expect_identical(extraction$extraction_mode, "BOTH_PRIMERS_SINGLE_READ")
    expect_true(extraction$full_insert_recovered)
    expect_identical(as.character(extraction$sequences_by_round$round_00), c("ACGTACGT", "TTTT"))
    expect_identical(extraction$input_reads, c(round_00 = 3))
    expect_identical(extraction$output_reads, c(round_00 = 2))
})

test_that("selexprep_extract supports one-sided modes and reverse orientation", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    five_prime <- selexprep_extract(
        list(round_00 = c(paste0(primer_5p, "ACGT"), "CCCC")),
        .report_with_mode("FIVE_PRIME_ONLY", primer_3p = NULL)
    )
    expect_identical(as.character(five_prime$sequences_by_round$round_00), "ACGT")

    reverse_input <- paste0(reverse_complement(primer_3p), "AACC", reverse_complement(primer_5p))
    reverse <- selexprep_extract(
        list(round_00 = reverse_input),
        .report_with_mode("BOTH_PRIMERS_SINGLE_READ", orientation = "REVERSE")
    )
    expect_identical(as.character(reverse$sequences_by_round$round_00), "GGTT")
})

test_that("mixed orientation records strand evidence without flipping reads", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    forward <- paste0(primer_5p, "ACGT", primer_3p)
    reverse <- reverse_complement(forward)
    report <- .report_with_mode("BOTH_PRIMERS_SINGLE_READ", orientation = "MIXED")

    extraction <- selexprep_extract(list(round_00 = c(forward, reverse)), report)

    expect_identical(as.character(extraction$sequences_by_round$round_00), "ACGT")
    expect_identical(extraction$strand_distribution$round_00,
        c(forward = 1L, reverse = 1L, ambiguous = 0L)
    )
})

test_that("selexprep_extract keeps paired split-primer sides separate", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    extraction <- selexprep_extract(
        list(round_00 = c(paste0(primer_5p, "ACGT"), "AAAA")),
        .report_with_mode("PAIRED_END_SPLIT_PRIMERS"),
        paired_mate_streams = list(round_00 = c(paste0(reverse_complement(primer_3p), "TGCA"), "CCCC"))
    )

    round <- extraction$sequences_by_round$round_00
    expect_identical(as.character(round$r1), "ACGT")
    expect_identical(as.character(round$r2), "TGCA")
    expect_false(extraction$full_insert_recovered)
})

test_that("selexprep_extract refuses unsafe reports", {
    unsafe <- selexprep:::.unable_library_report("R1", 42, "No sequences provided")
    expect_error(
        selexprep_extract(list(round_00 = "ACGT"), unsafe),
        "does not permit extraction"
    )
})
