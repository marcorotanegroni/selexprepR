.ena_fixture <- paste(
    paste(.ena_filereport_fields, collapse = "\t"),
    paste(
        "SRR000001", "PRJNA000001", "Example HT-SELEX", "SELEX", "TRANSCRIPTOMIC",
        "library_1", "experiment_1", "Example Round 1", "SRS000001",
        "1000", "80000", "aaa;bbb", "100;200",
        "ftp.sra.ebi.ac.uk/vol1/SRR000001_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/SRR000001_2.fastq.gz",
        sep = "\t"
    ),
    sep = "\n"
)

test_that("ENA parser preserves inspection metadata and builds a safe dry-run plan", {
    inspection <- .parse_ena_filereport(.ena_fixture, "PRJNA000001")
    result <- fetch_selexprep_reads(inspection, tempfile("selexprep-fetch-"), dry_run = TRUE)

    expect_s3_class(inspection, "selexprep_inspection")
    expect_identical(inspection$bioproject_id, "PRJNA000001")
    expect_identical(inspection$runs$fastq_urls[[1L]][[2L]],
        "ftp.sra.ebi.ac.uk/vol1/SRR000001_2.fastq.gz"
    )
    expect_s3_class(result, "selexprep_fetch_result")
    expect_true(result$dry_run)
    expect_true(all(result$plan$round_unambiguous))
    expect_true(all(result$plan$round_number == 1L))
    expect_identical(result$plan$url,
        c(
            "https://ftp.sra.ebi.ac.uk/vol1/SRR000001_1.fastq.gz",
            "https://ftp.sra.ebi.ac.uk/vol1/SRR000001_2.fastq.gz"
        )
    )
})

test_that("fetch refuses ambiguous round metadata unless acquisition-only is explicit", {
    inspection <- .parse_ena_filereport(.ena_fixture, "PRJNA000001")
    inspection$runs$sample_title <- "Example Round 1 Round 2"

    expect_error(
        fetch_selexprep_reads(inspection, tempfile("selexprep-fetch-"), dry_run = TRUE),
        "no unambiguous SELEX round assignment"
    )
    acquisition_only <- fetch_selexprep_reads(
        inspection,
        tempfile("selexprep-fetch-"),
        dry_run = TRUE,
        require_assigned_rounds = FALSE
    )
    expect_false(all(acquisition_only$plan$round_unambiguous))
})

test_that("fetch validates the download timeout before planning", {
    inspection <- .parse_ena_filereport(.ena_fixture, "PRJNA000001")

    expect_error(
        fetch_selexprep_reads(
            inspection,
            tempfile("selexprep-fetch-"),
            timeout_seconds = 0
        ),
        "positive number"
    )
})

test_that("ENA parser rejects missing fields and empty reports", {
    expect_error(.parse_ena_filereport("run_accession\nSRR000001", "PRJNA000001"), "missing required")
    expect_error(.parse_ena_filereport(paste(.ena_filereport_fields, collapse = "\t"), "PRJNA000001"),
        "no read-run records"
    )
})
