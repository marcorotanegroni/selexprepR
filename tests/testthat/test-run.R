.run_synthetic_pool <- function(round, primer_5p, primer_3p = NULL, random_length = 20L) {
    bases <- c("A", "C", "G", "T")
    vapply(0:499, function(index) {
        random_region <- paste0(
            bases[((index + round * 1000L) * 7L + seq.int(0L, random_length - 1L) * 13L) %% 4L + 1L],
            collapse = ""
        )
        paste0(primer_5p, random_region, if (is.null(primer_3p)) "" else primer_3p)
    }, character(1))
}

test_that("run_selexprep returns a self-describing SummarizedExperiment", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    pools <- stats::setNames(lapply(0:2, function(round) {
        .run_synthetic_pool(round, primer_5p, primer_3p)
    }), sprintf("round_%02d", 0:2))

    result <- run_selexprep(pools, accession = "PRJNA000001", low_total_reads = 0)

    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(
        unname(Matrix::colSums(SummarizedExperiment::assay(result, "counts"))),
        c(500, 500, 500)
    )
    expect_identical(S4Vectors::metadata(result)$library_report$extraction_mode,
        "BOTH_PRIMERS_SINGLE_READ"
    )
    expect_s3_class(S4Vectors::metadata(result)$extraction, "selexprep_extraction")
    expect_s3_class(S4Vectors::metadata(result)$qc, "selexprep_qc")
    expect_s3_class(S4Vectors::metadata(result)$manifest, "selexprep_manifest")
    expect_identical(S4Vectors::metadata(result)$manifest$accession, "PRJNA000001")
    expect_false(S4Vectors::metadata(result)$counting_skipped)
})

test_that("run_selexprep does not count half-inserts from a split pair", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    r1 <- stats::setNames(lapply(0:2, function(round) {
        .run_synthetic_pool(round, primer_5p, random_length = 80L)
    }), sprintf("round_%02d", 0:2))
    r2 <- stats::setNames(lapply(0:2, function(round) {
        .run_synthetic_pool(round, reverse_complement(primer_3p), random_length = 80L)
    }), sprintf("round_%02d", 0:2))

    result <- run_selexprep(r1, read_source = "R1_AND_R2", paired_mate_streams = r2)

    expect_identical(dim(result), c(0L, 3L))
    expect_true(S4Vectors::metadata(result)$counting_skipped)
    expect_null(S4Vectors::metadata(result)$qc)
    expect_identical(S4Vectors::metadata(result)$library_report$required_action,
        "READ_MERGING_RECOMMENDED"
    )
})

test_that("run_selexprep exposes extraction provenance to QC", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    forward <- paste0(primer_5p, "ACGT", primer_3p)
    reverse <- reverse_complement(forward)
    reads <- c(rep(forward, 400L), rep(reverse, 125L))

    experiment <- run_selexprep(
        list(round_00 = reads),
        low_total_reads = 0L
    )
    flags <- as.character(S4Vectors::metadata(experiment)$qc$flags$name)

    expect_true("strand_mix" %in% flags)
})

test_that("run_selexprep records manually reviewed primer overrides", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    reads <- rep(paste0(primer_5p, "ACGT", primer_3p), 100L)

    experiment <- run_selexprep(
        list(round_00 = reads),
        low_total_reads = 0L,
        primer_5p = primer_5p,
        primer_3p = primer_3p
    )
    metadata <- S4Vectors::metadata(experiment)

    expect_identical(metadata$library_report$status, "LOW")
    expect_identical(
        metadata$library_report$extraction_mode,
        "BOTH_PRIMERS_SINGLE_READ"
    )
    expect_identical(
        metadata$manifest$parameters[["manual_primer_override"]],
        "true"
    )
    expect_equal(
        unname(Matrix::colSums(SummarizedExperiment::assay(experiment))),
        100
    )
})
