test_that("low-yield diagnostics use both Python v0.4.6 thresholds", {
    low_yield <- selexprepR:::.low_yield_rounds
    expect_length(low_yield(c(empty = 0), c(empty = 0)), 0L)
    expect_length(low_yield(c(a = 100, b = 100), c(a = 100, b = 50)), 0L)
    expect_identical(names(low_yield(c(a = 100, b = 100), c(a = 100, b = 49))), "b")
    expect_identical(names(low_yield(c(a = 100, b = 100), c(a = 5, b = 5))), c("a", "b"))
    expect_match(unname(low_yield(c(a = 100), c(a = 5))), "5.0%; best round 5.0%", fixed = TRUE)
})

test_that("single and paired extraction warn and retain low-yield round evidence", {
    report <- read_library_report(test_path("fixtures", "library-report-both-primers.json"))
    good <- paste0(report$primer_5p, "ACGT", report$primer_3p)
    pools <- list(round_00 = rep(good, 10), round_01 = c(good, rep("CCCC", 9)))
    expect_warning(extraction <- selexprep_extract(pools, report), "primers may not fit")
    expect_identical(names(extraction$low_yield_rounds), "round_01")
    expect_identical(extraction$output_reads, c(round_00 = 10, round_01 = 1))
    expect_match(unname(extraction$low_yield_rounds), "kept 1 of 10 reads", fixed = TRUE)
    paired <- read_library_report(test_path("fixtures", "library-report-paired-split.json"))
    r1 <- list(round_00 = c(paste0(paired$primer_5p, "ACGT"), rep("AAAA", 3)))
    r2 <- list(round_00 = c(paste0(reverse_complement(paired$primer_3p), "TGCA"), rep("CCCC", 3)))
    expect_warning(split <- selexprep_extract(r1, paired, paired_mate_streams = r2), "round_00")
    expect_identical(names(split$low_yield_rounds), "round_00")
    expect_identical(as.character(split$sequences_by_round$round_00$r1), "ACGT")
    expect_false(split$full_insert_recovered)
})

test_that("lowercase RNA primers match while preserving insert case and alphabet", {
    report <- read_library_report(test_path("fixtures", "library-report-both-primers.json"))
    p5 <- tolower(chartr("T", "U", report$primer_5p))
    p3 <- tolower(chartr("T", "U", report$primer_3p))
    result <- selexprep_extract(list(round_00 = paste0(p5, "acgu", p3)), report)
    expect_identical(as.character(result$sequences_by_round$round_00), "acgu")
    expect_length(result$low_yield_rounds, 0L)
})

test_that("catalog values and source metadata are pinned to Python v0.4.6", {
    catalog <- selexprep_catalog()
    provenance <- S4Vectors::metadata(catalog)
    expect_identical(as.character(catalog$target_class[catalog$bioproject_id == "PRJEB114397"]), "small molecule")
    expect_identical(provenance$source_commit, "63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5")
    expect_identical(provenance$snapshot_date, "2026-10-02")
    expect_identical(digest::digest(system.file("extdata", "selexprep_catalog.csv", package = "selexprepR"),
        file = TRUE, algo = "sha256"), provenance$source_sha256)
})
