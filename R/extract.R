.reverse_complement_read <- function(sequence) {
    bases <- strsplit(sequence, "", fixed = TRUE)[[1L]]
    complements <- c(A = "T", C = "G", G = "C", T = "A", U = "A", N = "N", a = "t",
        c = "g", g = "c", t = "a", u = "a", n = "n")
    paste0(vapply(rev(bases), function(base) {
        complement <- complements[[base]]
        if (is.null(complement))
            "N" else complement
    }, character(1)), collapse = "")
}
.has_hamming_prefix <- function(sequence, primer) {
    sequence <- toupper(chartr("U", "T", sequence))
    primer <- toupper(chartr("U", "T", primer))
    nchar(sequence) >= nchar(primer) && .hamming_leq_one(substr(sequence, 1L, nchar(primer)),
        primer)
}
.normalise_adapter_alphabet <- function(sequence) {
    toupper(chartr("U", "T", sequence))
}
.adapter_match <- function(sequence, primer) {
    sequence <- .normalise_adapter_alphabet(sequence)
    primer <- .normalise_adapter_alphabet(primer)
    max_errors <- floor(nchar(primer) * 0.1)
    subject <- Biostrings::BString(sequence)
    pattern <- Biostrings::BString(primer)
    hits <- Biostrings::matchPattern(pattern, subject, max.mismatch = max_errors,
        with.indels = TRUE, fixed = TRUE)
    if (!length(hits)) {
        return(NULL)
    }
    starts <- Biostrings::start(hits)
    ends <- Biostrings::end(hits)
    errors <- Biostrings::neditAt(pattern, subject, at = starts, with.indels = TRUE,
        fixed = TRUE)
    widths <- ends - starts + 1L
    scores <- pmin(nchar(primer), widths) - 2L * errors
    best <- order(-scores, errors, starts, -ends, method = "radix")[[1L]]
    c(start = starts[[best]], end = ends[[best]])
}
.after_adapter_match <- function(sequence, match) {
    if (match[["end"]] >= nchar(sequence)) {
        return("")
    }
    substr(sequence, match[["end"]] + 1L, nchar(sequence))
}
.before_adapter_match <- function(sequence, match) {
    if (match[["start"]] <= 1L) {
        return("")
    }
    substr(sequence, 1L, match[["start"]] - 1L)
}
.trim_round <- function(sequences, mode, primer_5p, primer_3p) {
    kept <- switch(mode, BOTH_PRIMERS_SINGLE_READ = vapply(sequences, function(sequence) {
        match_5p <- .adapter_match(sequence, primer_5p)
        if (is.null(match_5p)) {
            return(NA_character_)
        }
        remainder <- .after_adapter_match(sequence, match_5p)
        match_3p <- .adapter_match(remainder, primer_3p)
        if (is.null(match_3p)) {
            return(NA_character_)
        }
        .before_adapter_match(remainder, match_3p)
    }, character(1)), FIVE_PRIME_ONLY = vapply(sequences, function(sequence) {
        match_5p <- .adapter_match(sequence, primer_5p)
        if (is.null(match_5p)) {
            return(NA_character_)
        }
        .after_adapter_match(sequence, match_5p)
    }, character(1)), THREE_PRIME_ONLY = vapply(sequences, function(sequence) {
        match_3p <- .adapter_match(sequence, primer_3p)
        if (is.null(match_3p)) {
            return(NA_character_)
        }
        .before_adapter_match(sequence, match_3p)
    }, character(1)), stop("Unsupported single-read extraction mode.", call. = FALSE))
    kept[!is.na(kept)]
}
.strand_distribution <- function(sequences, primer_5p, primer_3p) {
    reverse_3p <- if (is.null(primer_3p)) {
        NULL
    } else {
        reverse_complement(primer_3p)
    }
    forward <- reverse <- ambiguous <- 0L
    for (sequence in sequences) {
        if (!is.null(primer_5p) && .has_hamming_prefix(sequence, primer_5p)) {
            forward <- forward + 1L
        } else if (!is.null(reverse_3p) && .has_hamming_prefix(sequence, reverse_3p)) {
            reverse <- reverse + 1L
        } else {
            ambiguous <- ambiguous + 1L
        }
    }
    c(forward = forward, reverse = reverse, ambiguous = ambiguous)
}
.new_extraction <- function(mode, full_insert_recovered, sequences_by_round, input_reads,
    output_reads, strand_distribution = list()) {
    structure(list(extraction_mode = mode, full_insert_recovered = full_insert_recovered,
        sequences_by_round = sequences_by_round, input_reads = input_reads, output_reads = output_reads,
        strand_distribution = strand_distribution), class = c("selexprep_extraction",
        "list"))
}

.validate_primer_override <- function(primer, argument) {
    if (is.null(primer)) {
        return(NULL)
    }
    valid <- is.character(primer) && length(primer) == 1L &&
        !is.na(primer) && nzchar(primer)
    if (!valid || grepl("[^ACGTU]", toupper(primer))) {
        stop(
            sprintf(
                "`%s` must be one non-empty DNA or RNA sequence.",
                argument
            ),
            call. = FALSE
        )
    }
    toupper(primer)
}

.apply_primer_overrides <- function(report, primer_5p = NULL,
    primer_3p = NULL, extraction_mode = NULL) {
    requested <- !is.null(primer_5p) || !is.null(primer_3p) ||
        !is.null(extraction_mode)
    if (!requested) {
        return(report)
    }
    primer_5p <- .validate_primer_override(primer_5p, "primer_5p")
    primer_3p <- .validate_primer_override(primer_3p, "primer_3p")
    fields <- unclass(report)
    if (!is.null(primer_5p)) {
        fields$primer_5p <- primer_5p
        fields$variants_5p <- list()
    }
    if (!is.null(primer_3p)) {
        fields$primer_3p <- primer_3p
        fields$variants_3p <- list()
    }
    if (is.null(extraction_mode)) {
        extraction_mode <- fields$extraction_mode
        if (identical(extraction_mode, "UNABLE_TO_EXTRACT")) {
            extraction_mode <- if (!is.null(fields$primer_5p) &&
                !is.null(fields$primer_3p)) {
                "BOTH_PRIMERS_SINGLE_READ"
            } else if (!is.null(fields$primer_5p)) {
                "FIVE_PRIME_ONLY"
            } else if (!is.null(fields$primer_3p)) {
                "THREE_PRIME_ONLY"
            } else {
                "UNABLE_TO_EXTRACT"
            }
        }
    }
    valid_mode <- is.character(extraction_mode) &&
        length(extraction_mode) == 1L && !is.na(extraction_mode) &&
        extraction_mode %in% setdiff(.extraction_modes, "UNABLE_TO_EXTRACT")
    if (!valid_mode) {
        stop("`extraction_mode` must name an extractable mode.", call. = FALSE)
    }
    required <- switch(extraction_mode,
        BOTH_PRIMERS_SINGLE_READ = c("primer_5p", "primer_3p"),
        FIVE_PRIME_ONLY = "primer_5p",
        THREE_PRIME_ONLY = "primer_3p",
        PAIRED_END_SPLIT_PRIMERS = c("primer_5p", "primer_3p")
    )
    if (any(vapply(fields[required], is.null, logical(1)))) {
        stop(
            "The selected extraction mode requires additional primers.",
            call. = FALSE
        )
    }
    fields$extraction_mode <- extraction_mode
    fields$full_insert_recovered <- identical(
        extraction_mode,
        "BOTH_PRIMERS_SINGLE_READ"
    )
    fields$required_action <- if (identical(
        extraction_mode,
        "PAIRED_END_SPLIT_PRIMERS"
    )) {
        "READ_MERGING_RECOMMENDED"
    } else {
        "NONE"
    }
    fields$status <- "LOW"
    fields$confidence <- 0
    fields["failure_reason"] <- list(NULL)
    .new_library_report(fields)
}

#' Extract random regions from a SELEX library
#'
#' Uses a `selexprep_library_report` to remove the inferred constant regions.
#' The implementation operates directly on R sequence objects and preserves the
#' explicit safe-failure behavior of the original pipeline: reports marked as
#' unable to extract stop with an error unless the caller supplies a corrected
#' report. In split-primer paired-end mode, the two trimmed sides remain
#' separate; they are never silently merged.
#'
#' Full-primer matching follows the cutadapt defaults used by the Python
#' pipeline: adapters may occur away from the read boundary, substitutions
#' and indels are allowed up to a 10 percent error rate, and RNA `U` is
#' matched as DNA `T`. The returned insert retains the alphabet and case of
#' the input read.
#'
#' @param sequences_by_round A named list of character vectors or
#'   `Biostrings::XStringSet` objects containing R1 sequences.
#' @param library_report A `selexprep_library_report` from `selexprep_detect()`
#'   or `read_library_report()`.
#' @param paired_mate_streams Optional named R2 sequence pools, required for
#'   `PAIRED_END_SPLIT_PRIMERS` reports.
#' @param primer_5p Optional manually reviewed 5-prime primer override.
#' @param primer_3p Optional manually reviewed 3-prime primer override.
#' @param extraction_mode Optional manually reviewed extraction-mode override.
#'
#' @return A `selexprep_extraction`. Full and single-primer modes store one
#'   `Biostrings::BStringSet` per round; paired split mode stores an `r1` and
#'   `r2` `BStringSet` for every round.
#' @export
#'
#' @examples
#' primer_5p <- 'GGTAATACGACTCACTATAGGG'
#' primer_3p <- 'CCATGCATGCATGCATGCAT'
#' bases <- c('A', 'C', 'G', 'T')
#' random_regions <- vapply(0:499, function(i) {
#'     paste0(bases[(i * 7 + 0:15 * 13) %% 4 + 1], collapse = '')
#' }, character(1))
#' reads <- paste0(primer_5p, random_regions, primer_3p)
#' report <- selexprep_detect(list(round_00 = reads))
#' selexprep_extract(list(round_00 = reads[1:3]), report)
selexprep_extract <- function(sequences_by_round, library_report,
    paired_mate_streams = NULL, primer_5p = NULL, primer_3p = NULL,
    extraction_mode = NULL) {
    if (!inherits(library_report, "selexprep_library_report")) {
        stop("`library_report` must be a selexprep_library_report.", call. = FALSE)
    }
    library_report <- .apply_primer_overrides(
        library_report,
        primer_5p = primer_5p,
        primer_3p = primer_3p,
        extraction_mode = extraction_mode
    )
    if (identical(library_report$status, "UNABLE_TO_INFER") || identical(library_report$extraction_mode,
        "UNABLE_TO_EXTRACT")) {
        stop(paste("LibraryReport does not permit extraction;", "provide explicit corrected primers."),
            call. = FALSE)
    }
    pools <- .as_round_pools(sequences_by_round, "sequences_by_round")
    mode <- library_report$extraction_mode
    strand_distribution <- stats::setNames(lapply(pools, .strand_distribution,
        primer_5p = library_report$primer_5p, primer_3p = library_report$primer_3p),
        names(pools))
    if (identical(library_report$orientation, "REVERSE")) {
        pools <- lapply(pools, function(pool) {
            vapply(pool, .reverse_complement_read, character(1))
        })
    }
    if (identical(mode, "PAIRED_END_SPLIT_PRIMERS")) {
        if (is.null(paired_mate_streams)) {
            stop(paste("PAIRED_END_SPLIT_PRIMERS extraction requires", "`paired_mate_streams`."),
                call. = FALSE)
        }
        mate_pools <- .as_round_pools(paired_mate_streams, "paired_mate_streams")
        if (!identical(names(pools), names(mate_pools))) {
            stop(paste("`paired_mate_streams` must have the same round names as",
                "`sequences_by_round`."), call. = FALSE)
        }
        if (is.null(library_report$primer_5p) || is.null(library_report$primer_3p)) {
            stop("Split-primer extraction requires both primer sequences.", call. = FALSE)
        }
        primer_r2 <- reverse_complement(library_report$primer_3p)
        extracted <- lapply(names(pools), function(round_name) {
            r1 <- pools[[round_name]]
            r2 <- mate_pools[[round_name]]
            if (length(r1) != length(r2)) {
                stop(sprintf("Paired streams have different read counts in %s.",
                  round_name), call. = FALSE)
            }
            trimmed_r1 <- rep(NA_character_, length(r1))
            trimmed_r2 <- rep(NA_character_, length(r2))
            for (index in seq_along(r1)) {
                match_r1 <- .adapter_match(r1[[index]], library_report$primer_5p)
                match_r2 <- .adapter_match(r2[[index]], primer_r2)
                if (!is.null(match_r1) && !is.null(match_r2)) {
                  trimmed_r1[[index]] <- .after_adapter_match(r1[[index]], match_r1)
                  trimmed_r2[[index]] <- .after_adapter_match(r2[[index]], match_r2)
                }
            }
            keep <- !is.na(trimmed_r1)
            list(r1 = Biostrings::BStringSet(unname(trimmed_r1[keep])), r2 = Biostrings::BStringSet(unname(trimmed_r2[keep])))
        })
        names(extracted) <- names(pools)
        input_reads <- vapply(pools, length, numeric(1))
        output_reads <- vapply(extracted, function(round) length(round$r1), numeric(1))
        return(.new_extraction(mode, FALSE, extracted, input_reads, output_reads,
            strand_distribution))
    }
    if (identical(mode, "BOTH_PRIMERS_SINGLE_READ") && (is.null(library_report$primer_5p) ||
        is.null(library_report$primer_3p))) {
        stop("Both-primer extraction requires both primer sequences.", call. = FALSE)
    }
    if (identical(mode, "FIVE_PRIME_ONLY") && is.null(library_report$primer_5p)) {
        stop("Five-prime extraction requires a 5' primer.", call. = FALSE)
    }
    if (identical(mode, "THREE_PRIME_ONLY") && is.null(library_report$primer_3p)) {
        stop("Three-prime extraction requires a 3' primer.", call. = FALSE)
    }
    extracted <- lapply(pools, .trim_round, mode = mode, primer_5p = library_report$primer_5p,
        primer_3p = library_report$primer_3p)
    extracted <- lapply(extracted, function(sequences) {
        Biostrings::BStringSet(unname(sequences))
    })
    .new_extraction(mode, library_report$full_insert_recovered, extracted, vapply(pools,
        length, numeric(1)), vapply(extracted, length, numeric(1)), strand_distribution)
}
