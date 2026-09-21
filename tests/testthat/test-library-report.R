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

test_that("a zero-length inferred insert is refused", {
    dominant <- "TGCTTGGACTACATATGGTTGAGGGTTGTATGGAATTCTCGGGTGCCAAGG"
    pools <- stats::setNames(rep(list(rep(dominant, 600L)), 3L),
        sprintf("round_%02d", 0:2))

    report <- selexprep_detect(pools)
    expect_identical(report$status, "UNABLE_TO_INFER")
    expect_identical(report$extraction_mode, "UNABLE_TO_EXTRACT")
    expect_identical(report$required_action, "MANUAL_PRIMERS_REQUIRED")
    expect_null(report$primer_5p)
    expect_null(report$primer_3p)
    expect_match(report$failure_reason, "0 nt", fixed = TRUE)
})

test_that("outer-edge noise can be rescued without changing the insert boundary", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CGTGGTTACAGTCAGAGGACAGATT"
    flank_bases <- strsplit(primer_3p, "", fixed = TRUE)[[1L]]
    alternate <- c(A = "C", C = "G", G = "T", T = "A")
    bases <- c("A", "C", "G", "T")
    reads <- vapply(seq_len(600L), function(index) {
        flank <- vapply(seq_along(flank_bases), function(position) {
            base <- flank_bases[[position]]
            if (position <= 20L ||
                ((index * 7919L + position * 1543L) %% 1000L) < 600L) {
                base
            } else {
                alternate[[base]]
            }
        }, character(1))
        random_region <- paste0(
            bases[((index * 7L + seq_len(30L) * 13L) %% 4L) + 1L],
            collapse = ""
        )
        paste0(primer_5p, random_region, paste0(flank, collapse = ""))
    }, character(1))

    report <- selexprep_detect(rep(list(reads), 3L))
    expect_identical(report$extraction_mode, "BOTH_PRIMERS_SINGLE_READ")
    expect_identical(report$primer_3p, primer_3p)
    expect_identical(report$n_length_mode, 30)
    expect_gt(report$match_rate_3p, 0.7)
})

test_that("outer-edge core ignores isolated inner dips and clean flanks", {
    core <- selexprepR:::.high_support_core
    flank <- "CGTGGTTACAGTCAGAGGACAGATT"
    supports <- c(.921, .932, .931, .772, .932, .978, .923, .920, .911,
        .906, .903, .902, .898, .899, .884, .875, .976, .909, .901,
        .612, .483, .886, .896, .897, .749)

    expect_identical(core(flank, supports, FALSE), "CGTGGTTACAGTCAGAGGA")
    expect_null(core(flank, rep(.99, nchar(flank)), FALSE))
    clean_outer <- supports
    clean_outer[20:25] <- .95
    expect_null(core(flank, clean_outer, FALSE))
    expect_null(core(flank, c(rep(.99, 24L), .1), TRUE))
    expect_null(core("ACGTACGTACGTACGT",
        c(rep(.2, 5L), rep(.99, 11L)), TRUE))
    expect_identical(core("GACTACACTGCACTGCGTTAGAG",
        c(.55, .62, .71, rep(.99, 20L)), TRUE),
        "TACACTGCACTGCGTTAGAG")
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
    output_bytes <- readBin(path_a, "raw", n = file.info(path_a)$size)
    expect_false(as.raw(0x0d) %in% output_bytes)
    expect_false(endsWith(readChar(path_a, nchars = file.info(path_a)$size), "\n\n"))
    expect_equal(unclass(read_library_report(path_a)), unclass(report))
})
