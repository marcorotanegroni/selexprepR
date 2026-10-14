# Match the separators used by Python v0.4.6. Lookarounds retain separators
# between repeated indicators so R continues to expose every conflicting round.
# RV01/RV02 are ArrayExpress sample labels, never round indicators.
.round_patterns <- list(
    round_word_digit = "(?<![^[:space:]_\\-./])[Rr]ound[[:space:]_\\-]*(\\d+)(?=[[:space:]_\\-./]|$)",
    cycle_word_digit = "(?<![^[:space:]_\\-./])[Cc]ycle[[:space:]_\\-]*(\\d+)(?=[[:space:]_\\-./]|$)",
    iteration_digit = "(?<![^[:space:]_\\-./])[Ii]teration[[:space:]_\\-]*(\\d+)(?=[[:space:]_\\-./]|$)",
    pool_digit = "(?<![^[:space:]_\\-./])[Pp]ool[[:space:]_\\-]*(\\d+)(?=[[:space:]_\\-./]|$)",
    r_digit_boundary = "(?<![^[:space:]_\\-./])[Rr](\\d+)(?=[[:space:]_\\-./]|$)",
    c_digit_boundary = "(?<![^[:space:]_\\-./])[Cc](\\d+)(?=[[:space:]_\\-./]|$)",
    digit_r_suffix = "(?<![^[:space:]_\\-./])(\\d+)[Rr](?=[[:space:]_\\-./]|$)",
    glued_word_r_digit = "[A-Z][A-Za-z]{2,}R(\\d+)(?=[[:space:]_\\-./]|$)",
    digit_cyc_suffix = "(?i)(?<![^[:space:]_\\-./])(\\d+)[[:space:]_\\-]*cyc(?:les?)?(?=[[:space:]_\\-./]|$)"
)
.round_key_confidence <- function(key) {
    key <- tolower(trimws(key))
    key <- sub("^characteristics?[[:space:]]*\\[(.*)\\]$", "\\1", key, perl = TRUE)
    key <- gsub("[[:space:]\\-]+", "_", trimws(key), perl = TRUE)
    words <- strsplit(key, "_", fixed = TRUE)[[1L]]
    words <- words[nzchar(words)]
    round_words <- c("round", "cycle", "iteration", "r")
    vocabulary <- c("selex", "selection", "enrichment", "of", "number", "num", "no", "pcr")
    others <- words[!words %in% round_words]
    if (sum(words %in% round_words) != 1L || !all(others %in% vocabulary)) {
        return(NA_character_)
    }
    if ("pcr" %in% others) "MEDIUM" else "HIGH"
}
.round_attribute_value <- function(value) {
    pattern <- "(?i)^(?:(?:round|cycle|r|c)[[:space:]_\\-]*)?(\\d+)$"
    value <- trimws(value)
    if (is.na(value) || !grepl(pattern, value, perl = TRUE)) {
        return(NA_integer_)
    }
    suppressWarnings(as.integer(sub(pattern, "\\1", value, perl = TRUE)))
}
.round_matches <- function(text) {
    if (is.null(text) || !length(text) || is.na(text) || !nzchar(trimws(text))) {
        return(list(patterns = character(), values = integer()))
    }
    found_patterns <- character()
    found_values <- integer()
    for (name in names(.round_patterns)) {
        pattern <- .round_patterns[[name]]
        matches <- regmatches(text, gregexpr(pattern, text, perl = TRUE))[[1L]]
        if (!length(matches) || identical(matches, "")) {
            next
        }
        values <- suppressWarnings(as.integer(sub(pattern, "\\1", matches, perl = TRUE)))
        keep <- !is.na(values)
        found_patterns <- c(found_patterns, rep(name, sum(keep)))
        found_values <- c(found_values, values[keep])
    }
    list(patterns = found_patterns, values = found_values)
}
.as_round_attributes <- function(attributes) {
    if (is.null(attributes) || !length(attributes)) {
        return(character())
    }
    if (is.list(attributes)) {
        attributes <- unlist(attributes, use.names = TRUE)
    }
    if (is.null(names(attributes))) {
        return(character())
    }
    stats::setNames(as.character(attributes), names(attributes))
}
.target_hint <- function(sample_title) {
    if (is.null(sample_title) || !length(sample_title) || is.na(sample_title)) {
        return(NA_character_)
    }
    match <- regexec("^([A-Z][A-Za-z0-9_-]+)[[:space:]_-]+(?:[Rr]ound|[Rr]\\d|[Cc]ycle)",
        trimws(sample_title), perl = TRUE)
    found <- regmatches(trimws(sample_title), match)[[1L]]
    if (length(found) >= 2L)
        found[[2L]] else NA_character_
}
.text_round_assignment <- function(sample_title, library_name, experiment_title,
    design_description) {
    levels <- list(list(name = "sample_title", text = sample_title, confidence = "HIGH"),
        list(name = "library_name", text = library_name, confidence = "MEDIUM"),
        list(name = "experiment_title", text = experiment_title, confidence = "MEDIUM"),
        list(name = "design_description", text = design_description, confidence = "MEDIUM"))
    for (level in levels) {
        found <- .round_matches(level$text)
        if (!length(found$values)) {
            next
        }
        candidates <- unique(found$values)
        patterns <- unique(found$patterns)
        confidence <- if (length(candidates) == 1L)
            level$confidence else "MEDIUM"
        notes <- if (length(candidates) == 1L) {
            sprintf("matched '%s' in %s: '%s'", patterns[[1L]], level$name, substr(level$text,
                1L, 120L))
        } else {
            sprintf("conflicting candidates %s from patterns %s in %s: '%s'", paste(candidates,
                collapse = ","), paste(patterns, collapse = ","), level$name, substr(level$text,
                1L, 120L))
        }
        return(list(round_number = as.integer(candidates[[1L]]), confidence = confidence,
            source_field = level$name, matched_pattern = patterns[[1L]], candidates = as.integer(candidates),
            parser_notes = notes, target_hint = .target_hint(sample_title)))
    }
    NULL
}
.single_round_assignment <- function(run_accession, sample_title = "", library_name = "",
    experiment_title = "", design_description = "", sample_attributes = NULL, overrides = numeric()) {
    override <- if (run_accession %in% names(overrides))
        overrides[[run_accession]] else NA_integer_
    if (!is.null(override) && !is.na(override)) {
        return(list(round_number = as.integer(override), confidence = "HIGH", source_field = "seed_override",
            matched_pattern = "manual", candidates = as.integer(override), parser_notes = sprintf("overridden by supplied mapping: round=%d",
                as.integer(override)), target_hint = .target_hint(sample_title)))
    }
    attributes <- .as_round_attributes(sample_attributes)
    strong <- weak <- integer()
    for (index in seq_along(attributes)) {
        key <- names(attributes)[[index]]
        confidence <- .round_key_confidence(key)
        round <- .round_attribute_value(attributes[[index]])
        if (is.na(confidence) || is.na(round)) next
        if (confidence == "HIGH") {
            strong <- c(strong, stats::setNames(round, key))
        } else {
            weak <- c(weak, stats::setNames(round, key))
        }
    }
    text_record <- .text_round_assignment(sample_title, library_name,
        experiment_title, design_description)
    strong_rounds <- unique(unname(strong))
    weak_rounds <- unique(unname(weak))
    listed <- paste(sprintf("%s='%d'", names(c(strong, weak)), c(strong, weak)),
        collapse = "; ")
    if (length(strong_rounds) > 1L) {
        return(list(round_number = strong_rounds[[1L]], confidence = "MEDIUM",
            source_field = "sample_attributes", matched_pattern = "structured_attr",
            candidates = strong_rounds, parser_notes = sprintf("conflicting round attributes [%s]: %s",
                paste(strong_rounds, collapse = ", "), listed), target_hint = .target_hint(sample_title)))
    }
    if (length(strong_rounds) || (length(weak_rounds) && is.null(text_record))) {
        rounds <- if (length(strong_rounds)) strong_rounds else weak_rounds
        notes <- paste("structured attribute", listed)
        if (!length(strong_rounds)) {
            notes <- paste0(notes, "; the key can also mean amplification cycles, read as the round")
        } else if (length(weak_rounds) && !identical(weak_rounds, rounds)) {
            notes <- paste0(notes, "; ambiguous attribute reads as round [",
                paste(weak_rounds, collapse = ", "), "], not used")
        }
        fields <- list(sample_title = sample_title, library_name = library_name,
            experiment_title = experiment_title, design_description = design_description)
        for (name in names(fields)) {
            other <- unique(.round_matches(fields[[name]])$values)
            if (length(other) && !identical(other, rounds[1L])) {
                notes <- paste0(notes, "; ", name, " reads as round [",
                    paste(other, collapse = ", "), "] ('", substr(fields[[name]], 1L, 60L),
                    "'), the attribute is kept")
            }
        }
        return(list(round_number = rounds[[1L]], confidence = if (length(strong_rounds)) "HIGH" else "MEDIUM",
            source_field = "sample_attributes", matched_pattern = "structured_attr",
            candidates = rounds, parser_notes = notes, target_hint = .target_hint(sample_title)))
    }
    if (!is.null(text_record)) {
        if (length(weak_rounds) && !identical(weak_rounds, text_record$round_number)) {
            text_record$parser_notes <- paste0(text_record$parser_notes,
                "; ambiguous attribute ", listed, " not used over the text field")
        }
        return(text_record)
    }
    list(round_number = NA_integer_, confidence = "NONE", source_field = "none",
        matched_pattern = "none", candidates = integer(), parser_notes = "no round indicator found in supplied metadata fields",
        target_hint = .target_hint(sample_title))
}
.round_column <- function(runs, name, default = "") {
    if (!name %in% colnames(runs)) {
        return(rep(default, nrow(runs)))
    }
    runs[[name]]
}
#' Infer SELEX round assignments from public metadata
#'
#' Applies a conservative, ordered metadata cascade: structured sample
#' attributes, sample title, library name, experiment title, and design
#' description. A run remains unsafe when no round is found or when one field
#' contains conflicting candidates. The function records evidence rather than
#' using enrichment signals to repair metadata.
#'
#' Strong round attributes (including `Characteristics[selex cycle]`) take
#' precedence over text. Conflicting strong attributes remain ambiguous;
#' weak `pcr cycle` evidence never overrides an explicit text round.
#' ArrayExpress `RVxx` labels are sample numbers rather than round numbers.
#' Repeated conflicting indicators are all retained, including adjacent
#' `Round 1 Round 2` labels.
#'
#' @param runs A data frame or `S4Vectors::DataFrame` containing
#'   `run_accession` and optional metadata columns `sample_title`,
#'   `library_name`, `experiment_title`, `design_description`, and
#'   `sample_attributes`.
#' @param overrides Optional named integer vector mapping run accessions to
#'   manually curated round numbers.
#'
#' @return A `S4Vectors::DataFrame` with one conservative round assignment per
#'   run. `round_candidates` is a list-column and preserves ambiguity.
#' @export
#' @examples
#' runs <- data.frame(
#'     run_accession = c('SRR1', 'SRR2'),
#'     sample_title = c('Thrombin Round 1', 'Thrombin Round 2')
#' )
#' infer_selexprep_rounds(runs)
infer_selexprep_rounds <- function(runs, overrides = numeric()) {
    if (!inherits(runs, "DataFrame") && !is.data.frame(runs)) {
        stop("`runs` must be a DataFrame or data.frame.", call. = FALSE)
    }
    if (!"run_accession" %in% colnames(runs)) {
        stop("`runs` must contain a `run_accession` column.", call. = FALSE)
    }
    if (is.list(overrides)) {
        overrides <- unlist(overrides, use.names = TRUE)
    }
    if (!is.numeric(overrides) || (length(overrides) && (is.null(names(overrides)) ||
        any(!nzchar(names(overrides))) || anyNA(overrides) || any(overrides < 0) ||
        any(overrides != floor(overrides))))) {
        stop("`overrides` must be a named vector of non-negative integers.", call. = FALSE)
    }
    overrides <- stats::setNames(as.integer(overrides), names(overrides))
    accession <- as.character(runs$run_accession)
    if (anyNA(accession) || any(!nzchar(accession)) || anyDuplicated(accession)) {
        stop("`run_accession` must contain unique, non-missing strings.", call. = FALSE)
    }
    assignments <- lapply(seq_len(nrow(runs)), function(index) {
        .single_round_assignment(accession[[index]], sample_title = .round_column(runs,
            "sample_title")[[index]], library_name = .round_column(runs, "library_name")[[index]],
            experiment_title = .round_column(runs, "experiment_title")[[index]],
            design_description = .round_column(runs, "design_description")[[index]],
            sample_attributes = if ("sample_attributes" %in% colnames(runs)) {
                runs$sample_attributes[[index]]
            } else {
                NULL
            }, overrides = overrides)
    })
    S4Vectors::DataFrame(run_accession = accession, round_number = vapply(assignments,
        `[[`, integer(1), "round_number"), confidence = vapply(assignments, `[[`,
        character(1), "confidence"), source_field = vapply(assignments, `[[`, character(1),
        "source_field"), matched_pattern = vapply(assignments, `[[`, character(1),
        "matched_pattern"), round_candidates = I(lapply(assignments, `[[`, "candidates")),
        parser_notes = vapply(assignments, `[[`, character(1), "parser_notes"),
        target_hint = vapply(assignments, `[[`, character(1), "target_hint"))
}
.safe_round_assignments <- function(assignments) {
    !is.na(assignments$round_number) & lengths(assignments$round_candidates) ==
        1L
}
