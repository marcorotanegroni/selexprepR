.reverse_complement_read <- function(sequence) {
    bases <- strsplit(sequence, "", fixed = TRUE)[[1L]]
    complements <- c(A = "T", C = "G", G = "C", T = "A", U = "A", N = "N",
        a = "t", c = "g", g = "c", t = "a", u = "a", n = "n")
    paste0(vapply(rev(bases), function(base) {
        complement <- complements[[base]]
        if (is.null(complement)) "N" else complement
    }, character(1)),
        collapse = "")
}

.has_hamming_prefix <- function(sequence, primer) {
    nchar(sequence) >= nchar(primer) &&
        .hamming_leq_one(substr(sequence, 1L, nchar(primer)), primer)
}

.has_hamming_suffix <- function(sequence, primer) {
    sequence_length <- nchar(sequence)
    sequence_length >= nchar(primer) && .hamming_leq_one(
        substr(sequence, sequence_length - nchar(primer) + 1L, sequence_length), primer
    )
}

.trim_round <- function(sequences, mode, primer_5p, primer_3p) {
    kept <- switch(mode,
        BOTH_PRIMERS_SINGLE_READ = vapply(sequences, function(sequence) {
            if (!.has_hamming_prefix(sequence, primer_5p) || !.has_hamming_suffix(sequence, primer_3p)) {
                return(NA_character_)
            }
            substr(sequence, nchar(primer_5p) + 1L, nchar(sequence) - nchar(primer_3p))
        }, character(1)),
        FIVE_PRIME_ONLY = vapply(sequences, function(sequence) {
            if (!.has_hamming_prefix(sequence, primer_5p)) {
                return(NA_character_)
            }
            substr(sequence, nchar(primer_5p) + 1L, nchar(sequence))
        }, character(1)),
        THREE_PRIME_ONLY = vapply(sequences, function(sequence) {
            if (!.has_hamming_suffix(sequence, primer_3p)) {
                return(NA_character_)
            }
            substr(sequence, 1L, nchar(sequence) - nchar(primer_3p))
        }, character(1)),
        stop("Unsupported single-read extraction mode.", call. = FALSE)
    )
    kept[!is.na(kept)]
}

.strand_distribution <- function(sequences, primer_5p, primer_3p) {
    reverse_3p <- if (is.null(primer_3p)) NULL else reverse_complement(primer_3p)
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

.new_extraction <- function(mode, full_insert_recovered, sequences_by_round,
    input_reads, output_reads, strand_distribution = list()) {
    structure(list(
        extraction_mode = mode,
        full_insert_recovered = full_insert_recovered,
        sequences_by_round = sequences_by_round,
        input_reads = input_reads,
        output_reads = output_reads,
        strand_distribution = strand_distribution
    ), class = c("selexprep_extraction", "list"))
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
#' @param sequences_by_round A named list of character vectors or
#'   `Biostrings::XStringSet` objects containing R1 sequences.
#' @param library_report A `selexprep_library_report` from `selexprep_detect()`
#'   or `read_library_report()`.
#' @param paired_mate_streams Optional named R2 sequence pools, required for
#'   `PAIRED_END_SPLIT_PRIMERS` reports.
#'
#' @return A `selexprep_extraction`. Full and single-primer modes store one
#'   `Biostrings::BStringSet` per round; paired split mode stores an `r1` and
#'   `r2` `BStringSet` for every round.
#' @export
#'
#' @examples
#' primer_5p <- "GGTAATACGACTCACTATAGGG"
#' primer_3p <- "CCATGCATGCATGCATGCAT"
#' bases <- c("A", "C", "G", "T")
#' random_regions <- vapply(0:499, function(i) {
#'     paste0(bases[(i * 7 + 0:15 * 13) %% 4 + 1], collapse = "")
#' }, character(1))
#' reads <- paste0(primer_5p, random_regions, primer_3p)
#' report <- selexprep_detect(list(round_00 = reads))
#' selexprep_extract(list(round_00 = reads[1:3]), report)
selexprep_extract <- function(sequences_by_round, library_report, paired_mate_streams = NULL) {
    if (!inherits(library_report, "selexprep_library_report")) {
        stop("`library_report` must be a selexprep_library_report.", call. = FALSE)
    }
    if (identical(library_report$status, "UNABLE_TO_INFER") ||
        identical(library_report$extraction_mode, "UNABLE_TO_EXTRACT")) {
        stop("LibraryReport does not permit extraction; provide explicit corrected primers.",
            call. = FALSE)
    }
    pools <- .as_round_pools(sequences_by_round, "sequences_by_round")
    mode <- library_report$extraction_mode
    strand_distribution <- stats::setNames(lapply(pools, .strand_distribution,
        primer_5p = library_report$primer_5p,
        primer_3p = library_report$primer_3p
    ), names(pools))
    if (identical(library_report$orientation, "REVERSE")) {
        pools <- lapply(pools, function(pool) vapply(pool, .reverse_complement_read, character(1)))
    }

    if (identical(mode, "PAIRED_END_SPLIT_PRIMERS")) {
        if (is.null(paired_mate_streams)) {
            stop("PAIRED_END_SPLIT_PRIMERS extraction requires `paired_mate_streams`.", call. = FALSE)
        }
        mate_pools <- .as_round_pools(paired_mate_streams, "paired_mate_streams")
        if (!identical(names(pools), names(mate_pools))) {
            stop("`paired_mate_streams` must have the same round names as `sequences_by_round`.",
                call. = FALSE)
        }
        if (is.null(library_report$primer_5p) || is.null(library_report$primer_3p)) {
            stop("Split-primer extraction requires both primer sequences.", call. = FALSE)
        }
        primer_r2 <- reverse_complement(library_report$primer_3p)
        extracted <- lapply(names(pools), function(round_name) {
            r1 <- pools[[round_name]]
            r2 <- mate_pools[[round_name]]
            if (length(r1) != length(r2)) {
                stop(sprintf("Paired streams have different read counts in %s.", round_name), call. = FALSE)
            }
            keep <- vapply(seq_along(r1), function(index) {
                .has_hamming_prefix(r1[[index]], library_report$primer_5p) &&
                    .has_hamming_prefix(r2[[index]], primer_r2)
            }, logical(1))
            list(
                r1 = Biostrings::BStringSet(unname(substr(r1[keep], nchar(library_report$primer_5p) + 1L,
                    nchar(r1[keep])))),
                r2 = Biostrings::BStringSet(unname(substr(r2[keep], nchar(primer_r2) + 1L, nchar(r2[keep]))))
            )
        })
        names(extracted) <- names(pools)
        input_reads <- vapply(pools, length, numeric(1))
        output_reads <- vapply(extracted, function(round) length(round$r1), numeric(1))
        return(.new_extraction(mode, FALSE, extracted, input_reads, output_reads, strand_distribution))
    }

    if (identical(mode, "BOTH_PRIMERS_SINGLE_READ") &&
        (is.null(library_report$primer_5p) || is.null(library_report$primer_3p))) {
        stop("Both-primer extraction requires both primer sequences.", call. = FALSE)
    }
    if (identical(mode, "FIVE_PRIME_ONLY") && is.null(library_report$primer_5p)) {
        stop("Five-prime extraction requires a 5' primer.", call. = FALSE)
    }
    if (identical(mode, "THREE_PRIME_ONLY") && is.null(library_report$primer_3p)) {
        stop("Three-prime extraction requires a 3' primer.", call. = FALSE)
    }
    extracted <- lapply(pools, .trim_round, mode = mode,
        primer_5p = library_report$primer_5p, primer_3p = library_report$primer_3p)
    extracted <- lapply(extracted, function(sequences) Biostrings::BStringSet(unname(sequences)))
    .new_extraction(
        mode, library_report$full_insert_recovered, extracted,
        vapply(pools, length, numeric(1)), vapply(extracted, length, numeric(1)),
        strand_distribution
    )
}
