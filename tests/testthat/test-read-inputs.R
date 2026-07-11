.write_input_fastq <- function(path, sequences) {
    connection <- if (grepl("\\.gz$", path)) {
        gzfile(path, open = "wt")
    } else {
        file(path, open = "wt")
    }
    on.exit(close(connection), add = TRUE)
    records <- unlist(lapply(seq_along(sequences), function(index) {
        c(
            paste0("@read_", index), sequences[[index]], "+",
            paste(rep("I", nchar(sequences[[index]])), collapse = "")
        )
    }), use.names = FALSE)
    writeLines(records, connection)
    path
}

.input_synthetic_pool <- function(primer_5p, primer_3p) {
    bases <- c("A", "C", "G", "T")
    vapply(0:499, function(index) {
        random_region <- paste0(
            bases[(index * 7L + 0:15 * 13L) %% 4L + 1L],
            collapse = ""
        )
        paste0(primer_5p, random_region, primer_3p)
    }, character(1))
}

test_that("local FASTQ inputs are grouped by numeric round", {
    round_0 <- .write_input_fastq(
        tempfile(fileext = ".fastq"), c("ACGT", "TGCA")
    )
    round_1 <- .write_input_fastq(
        tempfile(fileext = ".fastq.gz"), c("AAAA", "CCCC", "GGGG")
    )

    inputs <- read_selexprep_inputs(
        c(round_00 = round_0, round_01 = round_1)
    )

    expect_s3_class(inputs, "selexprep_read_inputs")
    expect_identical(inputs$read_source, "R1")
    expect_null(inputs$paired_mate_streams)
    expect_identical(
        names(inputs$sequences_by_round), c("round_00", "round_01")
    )
    expect_identical(
        unname(as.character(inputs$sequences_by_round$round_01)),
        c("AAAA", "CCCC", "GGGG")
    )
    expect_s4_class(inputs$files, "DataFrame")
    expect_true(all(nchar(inputs$files$md5) == 32L))
})

test_that("completed fetch results preserve paired layout", {
    r1 <- .write_input_fastq(
        tempfile(pattern = "SRR000001_1_", fileext = ".fastq.gz"),
        c("AAAA", "CCCC")
    )
    r2 <- .write_input_fastq(
        tempfile(pattern = "SRR000001_2_", fileext = ".fastq.gz"),
        c("TTTT", "GGGG")
    )
    plan <- S4Vectors::DataFrame(
        run_accession = rep("SRR000001", 2L),
        round_number = rep(1L, 2L),
        round_unambiguous = rep(TRUE, 2L),
        file_name = c("SRR000001_1.fastq.gz", "SRR000001_2.fastq.gz")
    )
    fetch <- structure(list(
        plan = plan,
        downloaded_files = c(r1, r2),
        dry_run = FALSE
    ), class = c("selexprep_fetch_result", "list"))

    inputs <- read_selexprep_inputs(fetch)

    expect_identical(inputs$read_source, "R1_AND_R2")
    expect_identical(
        unname(as.character(inputs$sequences_by_round$round_01)),
        c("AAAA", "CCCC")
    )
    expect_identical(
        unname(as.character(inputs$paired_mate_streams$round_01)),
        c("TTTT", "GGGG")
    )
    expect_identical(
        as.character(inputs$files$run_accession),
        rep("SRR000001", 2L)
    )
})

test_that("input bridge rejects unsafe and incomplete layouts", {
    path <- .write_input_fastq(tempfile(fileext = ".fastq"), "ACGT")
    dry_run <- structure(list(
        plan = S4Vectors::DataFrame(),
        downloaded_files = character(),
        dry_run = TRUE
    ), class = c("selexprep_fetch_result", "list"))

    expect_error(read_selexprep_inputs(dry_run), "dry-run")
    expect_error(read_selexprep_inputs(path), "Supply `round`")
    expect_error(
        read_selexprep_inputs(c(path, path), round = c(0, 1)),
        "only once"
    )
    expect_error(
        read_selexprep_inputs(path, round = 0, mate = "R2"),
        "both R1 and R2"
    )
})

test_that("paired streams must contain equal numbers of reads", {
    r1 <- .write_input_fastq(
        tempfile(pattern = "sample_R1_", fileext = ".fastq"),
        c("AAAA", "CCCC")
    )
    r2 <- .write_input_fastq(
        tempfile(pattern = "sample_R2_", fileext = ".fastq"), "TTTT"
    )

    expect_error(
        read_selexprep_inputs(
            c(r1, r2), round = c(0, 0), mate = c("R1", "R2")
        ),
        "different read counts"
    )
})

test_that("FASTQ files run without manually assembled read lists", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    sequences <- .input_synthetic_pool(primer_5p, primer_3p)
    path <- .write_input_fastq(tempfile(fileext = ".fastq.gz"), sequences)

    experiment <- run_selexprep_files(
        path, round = 0, low_total_reads = 0
    )

    expect_s4_class(experiment, "SummarizedExperiment")
    expect_identical(colnames(experiment), "round_00")
    expect_s4_class(
        S4Vectors::metadata(experiment)$input_files, "DataFrame"
    )
    expect_identical(
        S4Vectors::metadata(experiment)$input_files$reads_loaded, 500L
    )
    manifest <- S4Vectors::metadata(experiment)$manifest
    expect_length(manifest$input_sha256, 1L)
    expect_match(names(manifest$input_sha256), "^round_00/R1/001_")
    expect_identical(
        unname(manifest$input_sha256),
        S4Vectors::metadata(experiment)$input_files$sha256
    )
})
