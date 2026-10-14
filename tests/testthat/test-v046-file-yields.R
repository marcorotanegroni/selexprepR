.file_yield_fastq <- function(sequences, pattern = "yield_") {
    path <- tempfile(pattern = pattern, fileext = ".fastq")
    records <- unlist(lapply(seq_along(sequences), function(index) {
        c(paste0("@read_", index), sequences[[index]], "+",
            strrep("I", nchar(sequences[[index]])))
    }), use.names = FALSE)
    writeLines(as.character(records), path)
    path
}

.file_yield_reads <- function(p5, p3 = "", width = 16L, n = 500L) {
    bases <- c("A", "C", "G", "T")
    vapply(seq_len(n), function(index) {
        random <- paste0(bases[(index * 7L + seq_len(width) * 13L) %% 4L + 1L],
            collapse = "")
        paste0(p5, random, p3)
    }, character(1))
}

test_that("file diagnostics reveal a bad FASTQ masked by the same round", {
    p5 <- "GGTAATACGACTCACTATAGGG"
    p3 <- "CCATGCATGCATGCATGCAT"
    bad <- .file_yield_fastq(rep("CCCC", 200L), "bad_")
    empty <- .file_yield_fastq(character(), "empty_")
    good <- .file_yield_fastq(.file_yield_reads(p5, p3, n = 600L), "good_")
    paths <- c(bad, empty, good)

    expect_warning(
        experiment <- run_selexprep_files(paths, round = 0L,
            primer_5p = p5, primer_3p = p3, low_total_reads = 0L),
        "input FASTQ.*kept 0 of 200 reads"
    )
    extraction <- S4Vectors::metadata(experiment)$extraction
    yields <- extraction$input_yields
    expect_s4_class(yields, "DataFrame")
    expect_identical(as.character(yields$path), normalizePath(paths, winslash = "/"))
    expect_identical(yields$input_reads, c(200L, 0L, 600L))
    expect_identical(yields$output_reads, c(0L, 0L, 600L))
    expect_identical(yields$fraction_kept, c(0, NA_real_, 1))
    expect_length(extraction$low_yield_rounds, 0L)
    expect_identical(names(extraction$low_yield_inputs), normalizePath(bad, winslash = "/"))
    expect_match(unname(extraction$low_yield_inputs), "best input 100.0%",
        fixed = TRUE)
    expect_null(attr(extraction, "kept_by_round", exact = TRUE))
    expect_equal(sum(SummarizedExperiment::assay(experiment)), 600)
    expect_length(S4Vectors::metadata(experiment)$manifest$input_sha256, 3L)
    expect_identical(S4Vectors::metadata(experiment)$input_files$reads_loaded,
        c(200L, 0L, 600L))
})

test_that("file yields use loaded reads after an explicit per-file cap", {
    p5 <- "GGTAATACGACTCACTATAGGG"
    p3 <- "CCATGCATGCATGCATGCAT"
    good <- .file_yield_fastq(.file_yield_reads(p5, p3, n = 600L), "good_")
    bad <- .file_yield_fastq(rep("CCCC", 600L), "bad_")
    expect_warning(
        experiment <- run_selexprep_files(c(good, bad), round = 0L,
            max_reads_per_file = 300L, primer_5p = p5, primer_3p = p3,
            low_total_reads = 0L),
        "input FASTQ"
    )
    yields <- S4Vectors::metadata(experiment)$extraction$input_yields
    expect_identical(yields$input_reads, c(300L, 300L))
    expect_identical(yields$output_reads, c(300L, 0L))
    expect_equal(sum(SummarizedExperiment::assay(experiment)), 300)
})

test_that("paired file yields count the jointly retained pairs once", {
    p5 <- "GGTAATACGACTCACTATAGGG"
    p3 <- "CCATGCATGCATGCATGCAT"
    bad_r1 <- .file_yield_fastq(.file_yield_reads(p5, width = 80L, n = 200L),
        "bad_R1_")
    bad_r2 <- .file_yield_fastq(rep("CCCC", 200L), "bad_R2_")
    good_r1 <- .file_yield_fastq(.file_yield_reads(p5, width = 80L, n = 600L),
        "good_R1_")
    good_r2 <- .file_yield_fastq(.file_yield_reads(reverse_complement(p3),
        width = 80L, n = 600L), "good_R2_")
    paths <- c(bad_r2, bad_r1, good_r2, good_r1)
    mates <- c("R2", "R1", "R2", "R1")
    expect_warning(
        experiment <- run_selexprep_files(paths, round = 0L, mate = mates,
            primer_5p = p5, primer_3p = p3,
            extraction_mode = "PAIRED_END_SPLIT_PRIMERS"),
        "1 input FASTQ.*kept 0 of 200 reads"
    )
    metadata <- S4Vectors::metadata(experiment)
    extraction <- metadata$extraction
    expect_identical(as.character(extraction$input_yields$path),
        normalizePath(c(bad_r1, good_r1), winslash = "/"))
    expect_identical(extraction$input_yields$input_reads, c(200L, 600L))
    expect_identical(extraction$input_yields$output_reads, c(0L, 600L))
    expect_identical(names(extraction$low_yield_inputs), normalizePath(bad_r1, winslash = "/"))
    expect_length(extraction$low_yield_rounds, 0L)
    expect_true(metadata$counting_skipped)
    expect_length(extraction$sequences_by_round$round_00$r1, 600L)
    expect_length(extraction$sequences_by_round$round_00$r2, 600L)
    expect_null(attr(extraction, "kept_by_round", exact = TRUE))
    expect_null(attr(extraction$sequences_by_round$round_00, "kept", exact = TRUE))
    expect_length(metadata$manifest$input_sha256, 4L)
})

test_that("single-read extraction excludes unused R2 files from yields", {
    p5 <- "GGTAATACGACTCACTATAGGG"
    p3 <- "CCATGCATGCATGCATGCAT"
    r1 <- .file_yield_fastq(.file_yield_reads(p5, p3), "both_R1_")
    r2 <- .file_yield_fastq(rep("CCCC", 500L), "unused_R2_")
    experiment <- run_selexprep_files(c(r1, r2), round = 0L,
        mate = c("R1", "R2"), primer_5p = p5, primer_3p = p3,
        extraction_mode = "BOTH_PRIMERS_SINGLE_READ", low_total_reads = 0L)
    metadata <- S4Vectors::metadata(experiment)
    expect_identical(as.character(metadata$extraction$input_yields$path),
        normalizePath(r1, winslash = "/"))
    expect_identical(metadata$extraction$input_yields$output_reads, 500L)
    expect_length(metadata$extraction$low_yield_inputs, 0L)
    expect_identical(nrow(metadata$input_files), 2L)
    expect_equal(sum(SummarizedExperiment::assay(experiment)), 500)
})

test_that("demultiplexed and direct pools retain no file-yield masks", {
    p5 <- "GGTAATACGACTCACTATAGGG"
    p3 <- "CCATGCATGCATGCATGCAT"
    reads <- .file_yield_reads(p5, p3)
    inputs <- selexprep_demultiplex(paste0("AAAAA", reads), c(AAAAA = 0L))
    experiment <- run_selexprep(inputs, low_total_reads = 0L)
    metadata <- S4Vectors::metadata(experiment)
    expect_null(metadata$extraction$input_yields)
    expect_null(metadata$extraction$low_yield_inputs)
    expect_null(attr(metadata$extraction, "kept_by_round", exact = TRUE))
    expect_identical(metadata$demultiplex$summary$assigned_reads, 500L)
    expect_equal(sum(SummarizedExperiment::assay(experiment)), 500)
    direct <- selexprep_extract(list(round_00 = reads), metadata$library_report)
    expect_null(attr(direct, "kept_by_round", exact = TRUE))
    expect_null(attr(direct$sequences_by_round$round_00, "kept", exact = TRUE))
})
