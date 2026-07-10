test_that("selexprep_count matches the Python count schema and values", {
    counts <- selexprep_count(c("ACGT", "ACGT", "ACGT", "GGGG", "TTTT"))

    expect_true(identical(colnames(counts), c("sequence", "reads", "rank", "rpm")))
    expect_true(identical(as.character(counts$sequence), c("ACGT", "GGGG", "TTTT")))
    expect_true(identical(counts$reads, c(3, 1, 1)))
    expect_true(identical(counts$rank, 1:3))
    expect_equal(counts$rpm, c(600000, 200000, 200000))
})

test_that("selexprep_count handles RNA and an empty round", {
    rna_counts <- selexprep_count(Biostrings::RNAStringSet(c("ACGU", "ACGU", "UUUU")))
    expect_true(identical(as.character(rna_counts$sequence), c("ACGU", "UUUU")))
    expect_true(identical(rna_counts$reads, c(2, 1)))

    empty <- selexprep_count(character())
    expect_identical(nrow(empty), 0L)
    expect_true(identical(colnames(empty), c("sequence", "reads", "rank", "rpm")))
})

test_that("selexprep_count validates its input", {
    expect_error(selexprep_count(1), "character vector")
    expect_error(selexprep_count(c("ACGT", NA_character_)), "missing values")
})

test_that("multi-round counts use a sparse SummarizedExperiment", {
    round_0 <- selexprep_count(c("AAAA", "AAAA", "CCCC"))
    round_1 <- selexprep_count(c("AAAA", "GGGG", "GGGG"))
    experiment <- selexprep:::.as_selexprep_experiment(
        list(round_00 = round_0, round_01 = round_1),
        accession = "PRJNA000001"
    )

    expect_s4_class(experiment, "SummarizedExperiment")
    expect_true(inherits(SummarizedExperiment::assay(experiment, "counts"), "dgCMatrix"))
    expect_true(identical(rownames(experiment), c("AAAA", "CCCC", "GGGG")))
    expect_true(identical(colnames(experiment), c("round_00", "round_01")))
    expect_equal(
        as.matrix(SummarizedExperiment::assay(experiment, "counts")),
        matrix(c(2, 1, 0, 1, 0, 2), nrow = 3, dimnames = list(
            c("AAAA", "CCCC", "GGGG"), c("round_00", "round_01")
        ))
    )
    expect_true(identical(as.character(SummarizedExperiment::rowData(experiment)$sequence),
        c("AAAA", "CCCC", "GGGG")
    ))
    expect_identical(S4Vectors::metadata(experiment)$accession, "PRJNA000001")
})
