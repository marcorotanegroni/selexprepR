# Expected round/confidence/source/candidates below were verified against
# selexprep v0.4.6 (63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5), parse_round.
.v046_round <- function(title = "", attributes = NULL, library = "") {
    runs <- S4Vectors::DataFrame(run_accession = "SRR1", sample_title = title,
        library_name = library, sample_attributes = I(list(attributes)))
    infer_selexprep_rounds(runs)
}

test_that("ArrayExpress RV sample numbers are not assigned as rounds", {
    for (fields in list(list(title = "RV04"), list(library = "RV04_s"))) {
        assignment <- do.call(.v046_round, fields)
        expect_identical(assignment$round_number, NA_integer_)
        expect_identical(assignment$confidence, "NONE")
        expect_false(.safe_round_assignments(assignment))
    }
    assignment <- .v046_round("RV04", c(construct = "ZFC1 of Sall4", "selex cycle" = "1"), "RV04_s")
    expect_identical(assignment$round_number, 1L)
    expect_identical(assignment$confidence, "HIGH")
    expect_identical(assignment$source_field, "sample_attributes")
})

test_that("structured round keys and values use the v0.4.6 vocabulary", {
    keys <- c("selex cycle", "Characteristics[selex cycle]", "SELEX-round",
        "selection cycle", "Selection cycle", "selex round", "round of selex enrichment",
        "round_selection", "round", "cycle number")
    for (key in keys) {
        assignment <- .v046_round(attributes = stats::setNames("6", key))
        expect_identical(assignment$round_number, 6L, info = key)
        expect_identical(assignment$source_field, "sample_attributes", info = key)
        expect_identical(assignment$confidence, "HIGH", info = key)
    }
    for (value in c("R3", "round 3", "Cycle-3", "3")) {
        expect_identical(.v046_round(attributes = c(round = value))$round_number, 3L, info = value)
    }
    for (key in c("cell cycle", "cycle threshold", "passage", "library_selection", "Selection", "round cycle")) {
        expect_identical(.v046_round(attributes = stats::setNames("3", key))$round_number,
            NA_integer_, info = key)
    }
})

test_that("all structured attributes are weighed independently of their order", {
    for (attributes in list(c("pcr cycle" = "25", "selex cycle" = "3"),
        c("selex cycle" = "3", "pcr cycle" = "25"))) {
        assignment <- .v046_round(attributes = attributes)
        expect_identical(assignment$round_number, 3L)
        expect_identical(assignment$confidence, "HIGH")
        expect_identical(assignment$round_candidates[[1L]], 3L)
        expect_match(assignment$parser_notes, "ambiguous attribute reads as round \\[25\\], not used")
    }
    for (attributes in list(c(round = "2", "selex cycle" = "3"),
        c("selex cycle" = "3", round = "2"))) {
        assignment <- .v046_round(attributes = attributes)
        expect_identical(sort(assignment$round_candidates[[1L]]), c(2L, 3L))
        expect_identical(assignment$confidence, "MEDIUM")
        expect_false(.safe_round_assignments(assignment))
        expect_match(assignment$parser_notes, "conflicting round attributes")
    }
    assignment <- .v046_round(attributes = c(round = "3", "selex cycle" = "Cycle3"))
    expect_identical(assignment$round_candidates[[1L]], 3L)
    expect_identical(assignment$confidence, "HIGH")
    expect_true(.safe_round_assignments(assignment))
})

test_that("strong attributes outrank text and PCR attributes do not", {
    strong <- .v046_round("Round 5", c("selex cycle" = "3"))
    expect_identical(strong$round_number, 3L)
    expect_identical(strong$source_field, "sample_attributes")
    expect_match(strong$parser_notes, "sample_title reads as round \\[5\\]")
    expect_false(grepl("reads as round", .v046_round("Round 3", c(round = "3"))$parser_notes))

    weak <- .v046_round("RV04", c("pcr cycle" = "1"))
    expect_identical(weak$round_number, 1L)
    expect_identical(weak$confidence, "MEDIUM")
    expect_true(.safe_round_assignments(weak))
    expect_match(weak$parser_notes, "amplification cycles")

    text <- .v046_round("Round 3", c("pcr cycle" = "25"))
    expect_identical(text$round_number, 3L)
    expect_identical(text$source_field, "sample_title")
    expect_match(text$parser_notes, "ambiguous attribute pcr cycle='25' not used")
    expect_false(grepl("not used", .v046_round("Round 3", c("pcr cycle" = "3"))$parser_notes))

    library <- .v046_round(attributes = c("pcr cycle" = "25"), library = "RAPT26-2R")
    expect_identical(library$round_number, 2L)
    expect_identical(library$source_field, "library_name")
})

test_that("text matching accepts the tagged Python separators and preserves ambiguity", {
    expect_identical(.v046_round(library = "SELEX_round26_N40")$round_number, 26L)
    expect_identical(.v046_round("DNAFOXR00")$round_number, 0L)
    expect_identical(.v046_round("Selex_7_cyc")$round_number, 7L)
    expect_identical(.v046_round("TFR1")$round_number, NA_integer_)
    expect_identical(.v046_round("(R3)")$round_number, NA_integer_)
    expect_identical(.v046_round("ROUND 3")$round_number, NA_integer_)
    conflict <- .v046_round("Round 3 / Cycle 7")
    expect_identical(conflict$round_candidates[[1L]], c(3L, 7L))
    expect_false(.safe_round_assignments(conflict))

    # R retains a deliberate conservative difference: Python consumes a shared
    # separator when matching, but R records every repeated conflicting label.
    repeated <- .v046_round("Round 1 Round 2")
    expect_identical(repeated$round_candidates[[1L]], c(1L, 2L))
    expect_false(.safe_round_assignments(repeated))
})
