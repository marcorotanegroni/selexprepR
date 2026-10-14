.v046_reference <- function() {
    jsonlite::fromJSON(test_path("fixtures", "v046-detection-reference.json"),
        simplifyVector = FALSE)
}
.v046_pools <- function(reference, descriptions) {
    pools <- lapply(descriptions, function(description) {
        inserts <- unlist(reference$insert_pools[[description$insert_pool]],
            use.names = FALSE)
        if (isTRUE(description$insert_edge_copies)) {
            copy <- (seq_along(inserts) - 1L) %% 5L < 3L
            inserts[copy] <- paste0(reference$core_5p,
                substr(inserts[copy], 13L, nchar(inserts[copy]) - 12L),
                reference$core_3p)
        }
        leading <- ""
        if (!is.null(description$leading_cycle)) {
            cycle <- strsplit(description$leading_cycle, "", fixed = TRUE)[[1L]]
            leading <- rep(cycle, length.out = length(inserts))
        }
        paste0(leading, description$prefix, inserts, description$suffix)
    })
    stats::setNames(pools, sprintf("round_%02d", seq_along(pools) - 1L))
}
.v046_expected_report <- function(payload) {
    as_named_numeric <- function(values) {
        flattened <- unlist(values, use.names = TRUE)
        result <- stats::setNames(as.numeric(flattened), names(flattened))
        if (length(result)) result[order(names(result), method = "radix")] else result
    }
    payload["known_adapter_hits"] <- list(as_named_numeric(payload$known_adapter_hits))
    payload["n_length_distribution"] <- list(as_named_numeric(payload$n_length_distribution))
    selexprepR:::.new_library_report(payload)
}

test_that("detection matches all fields of eight Python v0.4.6 reports", {
    reference <- .v046_reference()
    expect_identical(reference$source$tag, "v0.4.6")
    for (name in names(reference$cases)) {
        case <- reference$cases[[name]]
        pools <- .v046_pools(reference, case$rounds)
        mates <- if (is.null(case$mates)) NULL else .v046_pools(reference, case$mates)
        observed <- selexprep_detect(pools,
            read_source = if (is.null(mates)) "R1" else "R1_AND_R2",
            paired_mate_streams = mates)
        expected <- .v046_expected_report(case$report)
        expect_equal(unclass(observed), unclass(expected), tolerance = 1e-12,
            info = name)
    }
})

test_that("shared constants follow run-specific tag positions", {
    reference <- .v046_reference()
    pools <- .v046_pools(reference, reference$cases[["run-specific-tags"]]$rounds)
    constant <- selexprepR:::.shared_library_constant
    c5 <- "AGGACGATGCGGTACC"
    c3 <- "CAGACGACTCGCTGAGGATC"
    expect_identical(constant(paste0("ACTGAG", c5), pools, TRUE),
        list(constant = paste0("G", c5), offsets = c(5L, 5L, 6L)))
    expect_identical(constant(paste0(c3, "GCATG"), pools, FALSE),
        list(constant = c3, offsets = c(5L, 5L, 8L)))
})

test_that("unanchored rounds cannot trim the flank", {
    constant <- selexprepR:::.shared_library_constant
    flank <- "AGGACGATGCGGTACC"
    good <- rep(paste0(flank, "ACGTACGTACGTACGT"), 100L)
    absent <- rep(paste0(strrep("T", 16L), "ACGTACGTACGTACGT"), 100L)
    expect_identical(constant(flank, list(absent, good), TRUE),
        list(constant = flank, offsets = c(0L, 0L)))
    expect_identical(constant(flank, list(good, absent), TRUE),
        list(constant = flank, offsets = c(0L, NA_integer_)))
})

test_that("the inner core remains and half the reads can anchor a round", {
    constant <- selexprepR:::.shared_library_constant
    flank <- "AGGACGATGCGGTACC"
    core <- substr(flank, 5L, 16L)
    good <- paste0(flank, "ACGTACGTACGTACGT")
    changed <- rep(paste0("CCCC", core, "ACGTACGTACGTACGT"), 100L)
    absent <- paste0(strrep("T", 16L), "ACGTACGTACGTACGT")
    expect_identical(constant(flank,
        list(c(rep(good, 50L), rep(absent, 50L)), changed), TRUE),
        list(constant = core, offsets = c(4L, 4L)))
    expect_identical(constant(flank,
        list(c(rep(good, 49L), rep(absent, 51L)), changed), TRUE),
        list(constant = flank, offsets = c(0L, 0L)))
})

test_that("shared outer bases require the Python support and coverage floors", {
    constant <- selexprepR:::.shared_library_constant
    core <- "CGATGCGGTACC"
    flank <- paste0("A", core)
    a <- paste0("A", core, "ACGT")
    alternate <- paste0("C", core, "ACGT")
    bare <- paste0(core, "ACGT")
    expect_identical(constant(flank, list(c(rep(a, 11L), rep(alternate, 9L))), TRUE)$constant,
        flank)
    expect_identical(constant(flank, list(c(rep(a, 10L), rep(alternate, 10L))), TRUE)$constant,
        core)
    expect_identical(constant(flank, list(c(rep(a, 10L), rep(bare, 10L))), TRUE)$constant,
        flank)
    expect_identical(constant(flank, list(c(rep(a, 9L), rep(bare, 11L))), TRUE)$constant,
        core)
})

test_that("ambiguous core copies require eighty percent single-copy agreement", {
    constant <- selexprepR:::.shared_library_constant
    core <- "CGATGCGGTACC"
    flank <- paste0(strrep("A", 16L), core)
    good <- rep(paste0(flank, "ACGT"), 100L)
    inner <- paste0(flank, "ACGT")
    outer <- paste0(core, "ACGT")
    repeated <- paste0(core, "AAAA", core, "ACGT")
    at_threshold <- c(rep(inner, 80L), rep(outer, 20L), rep(repeated, 200L))
    below_threshold <- c(rep(inner, 79L), rep(outer, 21L), rep(repeated, 200L))
    expect_false(is.na(constant(flank, list(good, at_threshold), TRUE)$offsets[[2L]]))
    expect_identical(constant(flank, list(good, below_threshold), TRUE),
        list(constant = flank, offsets = c(0L, NA_integer_)))
})

test_that("a shortened constant that occurs twice in the flank is refused", {
    constant <- selexprepR:::.shared_library_constant
    core_5p <- "CGATGCGGTACC"
    flank_5p <- paste0(core_5p, "TTCA", core_5p)
    changed_5p <- paste0("A", substr(flank_5p, 2L, nchar(flank_5p)))
    first_5p <- paste0(c(rep(changed_5p, 70L), rep(flank_5p, 30L)), "ACGT")
    second_5p <- rep(paste0("CCCC", core_5p, "ACGT"), 100L)
    expect_identical(constant(flank_5p, list(first_5p, second_5p), TRUE),
        list(constant = flank_5p, offsets = c(0L, 0L)))

    core_3p <- "CAGACGACTCGC"
    flank_3p <- paste0(core_3p, "TTCA", core_3p)
    changed_3p <- paste0(substr(flank_3p, 1L, nchar(flank_3p) - 1L), "A")
    first_3p <- paste0("ACGT", c(rep(changed_3p, 70L), rep(flank_3p, 30L)))
    second_3p <- rep(paste0("ACGT", core_3p, "CCCC"), 100L)
    expect_identical(constant(flank_3p, list(first_3p, second_3p), FALSE),
        list(constant = flank_3p, offsets = c(0L, 0L)))
})

test_that("overlapping cores and flank search windows retain exact positions", {
    places <- selexprepR:::.core_places
    expect_identical(places("AAAA", "AA", TRUE, 4L), c(0L, 1L, 2L))
    expect_identical(places("AAAA", "AA", FALSE, 4L), c(2L, 1L, 0L))
    expect_identical(places("NNACGTNN", "ACGT", TRUE, 5L), integer())
    expect_identical(places("NNACGTNN", "ACGT", TRUE, 6L), 2L)
    expect_identical(places("NNACGTNN", "ACGT", FALSE, 5L), integer())
    expect_identical(places("NNACGTNN", "ACGT", FALSE, 6L), 2L)
    expect_identical(selexprepR:::.sequence_occurrences("AAAA", "AA"), 3L)
})

test_that("position and variant signals are measured at the constant offset", {
    primer <- "CGATGCGGTACC"
    prefix <- paste0(c("NN", "NNNNNNNN"), primer, "ACGT")
    suffix <- paste0("ACGT", primer, c("NN", "NNNNNNNN"))
    position <- selexprepR:::.position_consistency
    expect_identical(position(prefix, primer, TRUE, offset = 5L), 1)
    expect_identical(position(suffix, primer, FALSE, offset = 5L), 1)
    variants <- selexprepR:::.top_variants
    expect_identical(variants(prefix[[1L]], nchar(primer), TRUE, offset = 2L),
        list(list(primer, 1)))
    expect_identical(variants(suffix[[1L]], nchar(primer), FALSE, offset = 2L),
        list(list(primer, 1)))
})

test_that("lowercase RNA detection uses the same normalized pools as Python", {
    reference <- .v046_reference()
    case <- reference$cases[["shared-clean"]]
    pools <- .v046_pools(reference, case$rounds)
    rna <- lapply(pools, function(x) chartr("t", "u", tolower(x)))
    expect_identical(selexprepR:::.normalise_pool(c("uacgu", "UACGU")),
        c("TACGT", "TACGT"))
    expect_equal(unclass(selexprep_detect(rna)),
        unclass(.v046_expected_report(case$report)), tolerance = 1e-12)
})
