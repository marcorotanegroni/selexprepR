test_that("barcodes route, trim, and report single-end reads", {
    reads <- c(
        "AAAAACCCCC",
        "AAAATGGGGG",
        "TTTTTACGTA",
        "GGGGGAAAAA"
    )

    result <- selexprep_demultiplex(
        reads,
        c(AAAAA = 0L, TTTTT = 1L)
    )

    expect_s3_class(result, "selexprep_demultiplex")
    expect_identical(
        unname(as.character(result$sequences_by_round$round_00)),
        c("CCCCC", "GGGGG")
    )
    expect_identical(
        unname(as.character(result$sequences_by_round$round_01)),
        "ACGTA"
    )
    expect_identical(
        unname(as.character(result$unassigned$r1)),
        "GGGGGAAAAA"
    )
    expect_identical(result$summary$assigned_reads, c(2L, 1L))
})

test_that("paired mates remain synchronized and R2 is not trimmed", {
    r1 <- c("AAAAACCCCC", "TTTTTGGGGG", "CCCCCAAAAA")
    r2 <- c("R2_FIRST", "R2_SECOND", "R2_UNASSIGNED")

    result <- selexprep_demultiplex(
        r1,
        c(AAAAA = 0L, TTTTT = 1L),
        paired_mates = r2
    )

    expect_identical(result$read_source, "R1_AND_R2")
    expect_identical(
        unname(as.character(result$paired_mate_streams$round_00)),
        "R2_FIRST"
    )
    expect_identical(
        unname(as.character(result$paired_mate_streams$round_01)),
        "R2_SECOND"
    )
    expect_identical(
        unname(as.character(result$unassigned$r2)),
        "R2_UNASSIGNED"
    )
})

test_that("unsafe barcode codes and invalid arguments are rejected", {
    expect_error(
        selexprep_demultiplex(
            "AAAAACCCCC",
            c(AAAAA = 0L, AAAAT = 1L)
        ),
        "at least 3"
    )
    expect_error(
        selexprep_demultiplex(
            "AAAAACCCCC",
            c(AAAAA = 0L),
            paired_mates = character()
        ),
        "one mate"
    )
})

test_that("demultiplexed inputs run directly through the pipeline", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    reads <- rep(paste0("AAAAA", primer_5p, "ACGT", primer_3p), 500L)
    inputs <- selexprep_demultiplex(reads, c(AAAAA = 0L))

    experiment <- run_selexprep(inputs, low_total_reads = 0)

    expect_s4_class(experiment, "SummarizedExperiment")
    provenance <- S4Vectors::metadata(experiment)$demultiplex
    expect_identical(provenance$unassigned_reads, 0L)
    expect_identical(provenance$barcodes, c(AAAAA = 0L))
})
