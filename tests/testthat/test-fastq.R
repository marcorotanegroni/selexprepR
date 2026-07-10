test_that("FASTQ reader preserves RNA bases and supports gzip input", {
    path <- tempfile(fileext = ".fastq.gz")
    connection <- gzfile(path, open = "wt")
    writeLines(c("@read_1", "ACGU", "+", "IIII", "@read_2", "GGGG", "+", "JJJJ"), connection)
    close(connection)

    reads <- read_selexprep_fastq(path)

    expect_s4_class(reads, "BStringSet")
    expect_identical(unname(as.character(reads)), c("ACGU", "GGGG"))
    expect_identical(unname(as.character(read_selexprep_fastq(path, max_reads = 1))), "ACGU")
})

test_that("FASTQ reader validates the requested read cap", {
    path <- tempfile(fileext = ".fastq")
    writeLines(c("@read", "ACGT", "+", "IIII"), path)

    expect_error(read_selexprep_fastq(path, max_reads = 0), "positive integer")
    expect_error(read_selexprep_fastq("missing.fastq"), "existing FASTQ")
})
