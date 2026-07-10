test_that("selexprep count tables round-trip through portable formats", {
    counts <- selexprep_count(c("ACGT", "ACGT", "GGGG"))

    for (extension in c("tsv", "csv", "rds")) {
        path <- tempfile(fileext = paste0(".", extension))
        write_selexprep_counts(counts, path)
        restored <- read_selexprep_counts(path)
        expect_identical(as.character(restored$sequence), as.character(counts$sequence))
        expect_equal(restored$reads, counts$reads)
        expect_identical(restored$rank, counts$rank)
        expect_equal(restored$rpm, counts$rpm)
    }
})

test_that("selexprep count table reader rejects schema and RPM drift", {
    bad_schema <- tempfile(fileext = ".tsv")
    utils::write.table(
        data.frame(sequence = "ACGT", reads = 1, rank = 1, extra = 1),
        bad_schema,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    expect_error(read_selexprep_counts(bad_schema), "stable columns")

    bad_rpm <- tempfile(fileext = ".csv")
    utils::write.csv(
        data.frame(sequence = "ACGT", reads = 1, rank = 1, rpm = 1),
        bad_rpm,
        quote = FALSE,
        row.names = FALSE
    )
    expect_error(read_selexprep_counts(bad_rpm), "RPM")
})
