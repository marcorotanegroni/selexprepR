.synthetic_pool <- function(primer_5p = NULL, primer_3p = NULL, n = 500L,
    random_length = 30L, offset = 0L) {
    bases <- strsplit("ACGT", "", fixed = TRUE)[[1L]]
    vapply(seq.int(0L, n - 1L), function(index) {
        random_region <- paste0(
            bases[((index + offset) * 7L + seq.int(0L, random_length - 1L) * 13L) %% 4L + 1L],
            collapse = ""
        )
        paste0(
            if (is.null(primer_5p)) "" else primer_5p,
            random_region,
            if (is.null(primer_3p)) "" else primer_3p
        )
    }, character(1))
}

.three_round_pool <- function(primer_5p = NULL, primer_3p = NULL, n = 500L, random_length = 30L) {
    stats::setNames(lapply(0:2, function(round) {
        .synthetic_pool(primer_5p, primer_3p, n = n, random_length = random_length,
            offset = round * 1000L)
    }), sprintf("round_%02d", 0:2))
}

test_that("selexprep_detect matches the Python golden LibraryReport", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    observed <- selexprep_detect(.three_round_pool(primer_5p, primer_3p))
    reference <- read_library_report(test_path("fixtures", "library-report-both-primers.json"))

    expect_equal(unclass(observed), unclass(reference))
})

test_that("selexprep_detect matches Python golden reports for all classification modes", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    cases <- list(
        "five-prime-only" = list(
            observed = selexprep_detect(.three_round_pool(primer_5p, n = 1000L)),
            fixture = "library-report-five-prime-only.json"
        ),
        "three-prime-only" = list(
            observed = selexprep_detect(.three_round_pool(primer_3p = primer_3p, n = 1000L)),
            fixture = "library-report-three-prime-only.json"
        ),
        "paired-split" = list(
            observed = selexprep_detect(
                .three_round_pool(primer_5p, n = 1000L, random_length = 80L),
                read_source = "R1_AND_R2",
                paired_mate_streams = .three_round_pool(
                    reverse_complement(primer_3p), n = 1000L, random_length = 80L
                )
            ),
            fixture = "library-report-paired-split.json"
        ),
        "unable" = list(
            observed = selexprep_detect(list(round_00 = .synthetic_pool(n = 100L))),
            fixture = "library-report-unable.json"
        )
    )

    for (name in names(cases)) {
        reference <- read_library_report(test_path("fixtures", cases[[name]]$fixture))
        expect_equal(unclass(cases[[name]]$observed), unclass(reference), info = name)

        output <- tempfile(fileext = ".json")
        write_library_report(cases[[name]]$observed, output)
        reference_path <- test_path("fixtures", cases[[name]]$fixture)
        expect_identical(readBin(output, "raw", n = file.info(output)$size),
            readBin(reference_path, "raw", n = file.info(reference_path)$size),
            info = paste(name, "JSON"))
    }
})

test_that("paired split-primer libraries are surfaced without silent merging", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    r1 <- .three_round_pool(primer_5p, random_length = 80L)
    r2 <- .three_round_pool(reverse_complement(primer_3p), random_length = 80L)

    report <- selexprep_detect(r1, read_source = "R1_AND_R2", paired_mate_streams = r2)
    expect_identical(report$extraction_mode, "PAIRED_END_SPLIT_PRIMERS")
    expect_identical(report$required_action, "READ_MERGING_RECOMMENDED")
    expect_identical(report$primer_3p, primer_3p)
})

test_that("LibraryReport JSON is deterministic and round-trips", {
    report <- selexprep_detect(.three_round_pool(
        "GGTAATACGACTCACTATAGGG", "CCATGCATGCATGCATGCAT"
    ))
    path_a <- tempfile(fileext = ".json")
    path_b <- tempfile(fileext = ".json")

    hash_a <- write_library_report(report, path_a)
    hash_b <- write_library_report(report, path_b)
    expect_identical(hash_a, hash_b)
    expect_identical(readBin(path_a, "raw", n = file.info(path_a)$size),
        readBin(path_b, "raw", n = file.info(path_b)$size)
    )
    expect_false(endsWith(readChar(path_a, nchars = file.info(path_a)$size), "\n\n"))
    expect_equal(unclass(read_library_report(path_a)), unclass(report))
})
