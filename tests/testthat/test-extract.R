.report_with_mode <- function(mode, primer_5p = "GGTAATACGACTCACTATAGGG",
    primer_3p = "CCATGCATGCATGCATGCAT", orientation = "FORWARD") {
    fields <- unclass(read_library_report(test_path(
        "fixtures",
        "library-report-both-primers.json"
    )))
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
    selexprepR:::.new_library_report(fields)
}

test_that("selexprep_extract retains only fully linked-primer reads", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    good <- paste0(primer_5p, "ACGTACGT", primer_3p)
    one_mismatch <- paste0(
        "A",
        substr(primer_5p, 2L, nchar(primer_5p)),
        "TTTT",
        primer_3p
    )
    two_mismatch_primer <- primer_5p
    substr(two_mismatch_primer, 1L, 1L) <- "A"
    substr(two_mismatch_primer, 8L, 8L) <- "T"
    two_mismatches <- paste0(two_mismatch_primer, "CCCC", primer_3p)
    three_mismatch_primer <- two_mismatch_primer
    substr(three_mismatch_primer, 12L, 12L) <- "A"
    three_mismatches <- paste0(three_mismatch_primer, "AAAA", primer_3p)
    bad <- paste0(primer_5p, "GGGGGGGG")

    extraction <- selexprep_extract(
        list(round_00 = c(
            good,
            one_mismatch,
            two_mismatches,
            three_mismatches,
            bad
        )),
        .report_with_mode("BOTH_PRIMERS_SINGLE_READ")
    )

    expect_identical(extraction$extraction_mode, "BOTH_PRIMERS_SINGLE_READ")
    expect_true(extraction$full_insert_recovered)
    expect_identical(
        as.character(extraction$sequences_by_round$round_00),
        c("ACGTACGT", "TTTT", "CCCC")
    )
    expect_identical(extraction$input_reads, c(round_00 = 5))
    expect_identical(extraction$output_reads, c(round_00 = 3))
})

test_that("approximate linked matching handles indels away from boundaries", {
    primer_5p <- "ACGTCAGTACGATCGTACGA"
    primer_3p <- "TGCATCGACTGACGTAGCTA"
    deleted_5p <- paste0(
        substr(primer_5p, 1L, 6L),
        substr(primer_5p, 8L, nchar(primer_5p))
    )
    inserted_5p <- paste0(
        substr(primer_5p, 1L, 7L),
        "T",
        substr(primer_5p, 8L, nchar(primer_5p))
    )
    deleted_3p <- paste0(
        substr(primer_3p, 1L, 8L),
        substr(primer_3p, 10L, nchar(primer_3p))
    )
    inserted_3p <- paste0(
        substr(primer_3p, 1L, 9L),
        "A",
        substr(primer_3p, 10L, nchar(primer_3p))
    )
    reads <- c(
        paste0("NN", deleted_5p, "AAAA", primer_3p, "NN"),
        paste0("GG", inserted_5p, "CCCC", primer_3p, "TT"),
        paste0("TT", primer_5p, "GGGG", deleted_3p, "AA"),
        paste0("CC", primer_5p, "TTTT", inserted_3p, "GG")
    )

    extraction <- selexprep_extract(
        list(round_00 = reads),
        .report_with_mode(
            "BOTH_PRIMERS_SINGLE_READ",
            primer_5p,
            primer_3p
        )
    )

    expect_identical(
        as.character(extraction$sequences_by_round$round_00),
        c("AAAA", "CCCC", "GGGG", "TTTT")
    )
})

test_that("RNA primer notation matches DNA reads", {
    primer_5p_rna <- "ACGUUAGUACGAUCGUACGA"
    primer_3p_rna <- "UGCAUCGACUGACGUAGCUA"
    primer_5p_dna <- chartr("U", "T", primer_5p_rna)
    primer_3p_dna <- chartr("U", "T", primer_3p_rna)
    read <- paste0(primer_5p_dna, "ACGU", primer_3p_dna)

    extraction <- selexprep_extract(
        list(round_00 = read),
        .report_with_mode(
            "BOTH_PRIMERS_SINGLE_READ",
            primer_5p_rna,
            primer_3p_rna
        )
    )

    expect_identical(
        as.character(extraction$sequences_by_round$round_00),
        "ACGU"
    )
})

test_that("one-sided modes and reverse orientation are supported", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    five_prime <- selexprep_extract(
        list(round_00 = c(paste0(primer_5p, "ACGT"), "CCCC")),
        .report_with_mode("FIVE_PRIME_ONLY", primer_3p = NULL)
    )
    expect_identical(
        as.character(five_prime$sequences_by_round$round_00),
        "ACGT"
    )

    reverse_input <- paste0(
        reverse_complement(primer_3p),
        "AACC",
        reverse_complement(primer_5p)
    )
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
    report <- .report_with_mode(
        "BOTH_PRIMERS_SINGLE_READ",
        orientation = "MIXED"
    )

    extraction <- selexprep_extract(
        list(round_00 = c(forward, reverse)),
        report
    )

    expect_identical(
        as.character(extraction$sequences_by_round$round_00),
        "ACGT"
    )
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
        paired_mate_streams = list(round_00 = c(
            paste0(reverse_complement(primer_3p), "TGCA"),
            "CCCC"
        ))
    )

    round <- extraction$sequences_by_round$round_00
    expect_identical(as.character(round$r1), "ACGT")
    expect_identical(as.character(round$r2), "TGCA")
    expect_false(extraction$full_insert_recovered)
})

test_that("paired split extraction uses approximate primer coordinates", {
    primer_5p <- "ACGTCAGTACGATCGTACGA"
    primer_3p <- "TGCATCGACTGACGTAGCTA"
    primer_r2 <- reverse_complement(primer_3p)
    deleted_r1 <- paste0(
        substr(primer_5p, 1L, 6L),
        substr(primer_5p, 8L, nchar(primer_5p))
    )
    inserted_r2 <- paste0(
        substr(primer_r2, 1L, 9L),
        "A",
        substr(primer_r2, 10L, nchar(primer_r2))
    )

    extraction <- selexprep_extract(
        list(round_00 = paste0("NN", deleted_r1, "ACGT")),
        .report_with_mode(
            "PAIRED_END_SPLIT_PRIMERS",
            primer_5p,
            primer_3p
        ),
        paired_mate_streams = list(
            round_00 = paste0("GG", inserted_r2, "TGCA")
        )
    )

    round <- extraction$sequences_by_round$round_00
    expect_identical(as.character(round$r1), "ACGT")
    expect_identical(as.character(round$r2), "TGCA")
})

test_that("selexprep_extract refuses unsafe reports", {
    unsafe <- selexprepR:::.unable_library_report(
        "R1",
        42,
        "No sequences provided"
    )
    expect_error(
        selexprep_extract(list(round_00 = "ACGT"), unsafe),
        "does not permit extraction"
    )
})

test_that("reviewed primer overrides recover an unsafe report", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    unsafe <- selexprepR:::.unable_library_report(
        "R1",
        42,
        "No sequences provided"
    )

    extraction <- selexprep_extract(
        list(round_00 = paste0(primer_5p, "ACGT", primer_3p)),
        unsafe,
        primer_5p = primer_5p,
        primer_3p = primer_3p
    )

    expect_identical(
        as.character(extraction$sequences_by_round$round_00),
        "ACGT"
    )
    expect_identical(
        extraction$extraction_mode,
        "BOTH_PRIMERS_SINGLE_READ"
    )
    expect_identical(unsafe$extraction_mode, "UNABLE_TO_EXTRACT")
})
