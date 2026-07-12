.library_report_fields <- c("primer_5p", "primer_3p", "variants_5p", "variants_3p",
    "known_adapter_hits", "extraction_mode", "full_insert_recovered", "read_source",
    "required_action", "orientation", "n_length_mode", "n_length_distribution",
    "n_length_confidence", "match_rate_5p", "match_rate_3p", "position_consistency_5p",
    "position_consistency_3p", "read_fraction_used_for_inference", "sampling_seed",
    "confidence", "status", "failure_reason")
.extraction_modes <- c("BOTH_PRIMERS_SINGLE_READ", "FIVE_PRIME_ONLY", "THREE_PRIME_ONLY",
    "PAIRED_END_SPLIT_PRIMERS", "UNABLE_TO_EXTRACT")
.read_sources <- c("R1", "R2", "R1_AND_R2", "INTERLEAVED", "UNKNOWN")
.required_actions <- c("NONE", "MANUAL_PRIMERS_REQUIRED", "READ_MERGING_RECOMMENDED")
.orientations <- c("FORWARD", "REVERSE", "MIXED")
.statuses <- c("HIGH", "MEDIUM", "LOW", "UNABLE_TO_INFER")
.known_adapters <- c(TRUSEQ_R1 = "AGATCGGAAGAGC", NEXTERA = "CTGTCTCTTATACACATCT")
.adapter_probe_length <- 13L
.new_library_report <- function(fields) {
    if (!identical(names(fields), .library_report_fields)) {
        stop("A LibraryReport must contain the stable selexprep field set.", call. = FALSE)
    }
    for (name in c("primer_5p", "primer_3p", "failure_reason")) {
        value <- fields[[name]]
        if (!is.null(value) && (!is.character(value) || length(value) != 1L ||
            is.na(value))) {
            stop(sprintf("`%s` must be NULL or one non-missing string.", name),
                call. = FALSE)
        }
    }
    for (name in c("extraction_mode", "read_source", "required_action", "orientation",
        "status")) {
        allowed <- switch(name, extraction_mode = .extraction_modes, read_source = .read_sources,
            required_action = .required_actions, orientation = .orientations, status = .statuses)
        if (!is.character(fields[[name]]) || length(fields[[name]]) != 1L || !(fields[[name]] %in%
            allowed)) {
            stop(sprintf("`%s` contains an unsupported value.", name), call. = FALSE)
        }
    }
    if (!is.logical(fields$full_insert_recovered) || length(fields$full_insert_recovered) !=
        1L || is.na(fields$full_insert_recovered)) {
        stop("`full_insert_recovered` must be one non-missing logical value.",
            call. = FALSE)
    }
    numeric_fields <- c("n_length_confidence", "match_rate_5p", "match_rate_3p",
        "position_consistency_5p", "position_consistency_3p", "read_fraction_used_for_inference",
        "sampling_seed", "confidence")
    if (any(!vapply(fields[numeric_fields], function(x) is.numeric(x) && length(x) ==
        1L && !is.na(x), logical(1)))) {
        stop("LibraryReport numeric scalar fields must be non-missing scalars.",
            call. = FALSE)
    }
    if (!is.null(fields$n_length_mode) && (!is.numeric(fields$n_length_mode) ||
        length(fields$n_length_mode) != 1L)) {
        stop("`n_length_mode` must be NULL or a numeric scalar.", call. = FALSE)
    }
    if (!is.numeric(fields$n_length_distribution) || (length(fields$n_length_distribution) &&
        is.null(names(fields$n_length_distribution)))) {
        stop("`n_length_distribution` must be a named numeric vector.", call. = FALSE)
    }
    if (!is.numeric(fields$known_adapter_hits) || (length(fields$known_adapter_hits) &&
        is.null(names(fields$known_adapter_hits)))) {
        stop("`known_adapter_hits` must be a named numeric vector.", call. = FALSE)
    }
    if (!is.list(fields$variants_5p) || !is.list(fields$variants_3p)) {
        stop("LibraryReport variants must be lists.", call. = FALSE)
    }
    structure(fields, class = c("selexprep_library_report", "list"))
}
.unable_library_report <- function(read_source, sampling_seed, failure_reason,
    known_adapter_hits = numeric()) {
    .new_library_report(list(primer_5p = NULL, primer_3p = NULL, variants_5p = list(),
        variants_3p = list(), known_adapter_hits = known_adapter_hits, extraction_mode = "UNABLE_TO_EXTRACT",
        full_insert_recovered = FALSE, read_source = read_source, required_action = "MANUAL_PRIMERS_REQUIRED",
        orientation = "FORWARD", n_length_mode = NULL, n_length_distribution = numeric(),
        n_length_confidence = 0, match_rate_5p = 0, match_rate_3p = 0, position_consistency_5p = 0,
        position_consistency_3p = 0, read_fraction_used_for_inference = 0, sampling_seed = as.numeric(sampling_seed),
        confidence = 0, status = "UNABLE_TO_INFER", failure_reason = failure_reason))
}
.hamming_leq_one <- function(left, right) {
    if (nchar(left) != nchar(right)) {
        return(FALSE)
    }
    sum(strsplit(left, "", fixed = TRUE)[[1L]] != strsplit(right, "", fixed = TRUE)[[1L]]) <=
        1L
}
.first_mode <- function(values) {
    unique_values <- unique(values)
    counts <- tabulate(match(values, unique_values), nbins = length(unique_values))
    index <- which.max(counts)
    list(value = unique_values[[index]], count = counts[[index]])
}
.positional_consensus <- function(sequences, is_prefix, max_length, min_sequences) {
    consensus <- character()
    supports <- numeric()
    widths <- nchar(sequences)
    for (offset in seq.int(0L, max_length - 1L)) {
        keep <- widths > offset
        if (sum(keep) < min_sequences) {
            break
        }
        bases <- if (is_prefix) {
            substr(sequences[keep], offset + 1L, offset + 1L)
        } else {
            substr(sequences[keep], widths[keep] - offset, widths[keep] - offset)
        }
        mode <- .first_mode(bases)
        consensus <- c(consensus, mode$value)
        supports <- c(supports, mode$count/length(bases))
    }
    list(consensus = consensus, supports = supports)
}
.flank_boundary <- function(supports, min_length, confidence) {
    if (length(supports) < min_length) {
        return(0L)
    }
    boundary <- length(supports)
    for (i in seq.int(min_length + 1L, max(min_length + 1L, length(supports) -
        1L))) {
        if (i + 1L > length(supports)) {
            break
        }
        baseline <- stats::median(supports[(i - min_length):(i - 1L)])
        window <- supports[i:(i + 1L)]
        if (baseline >= 0.9 && supports[[i]] <= baseline - 0.12 && all(window <
            0.85)) {
            boundary <- i - 1L
            break
        }
    }
    floor <- min(confidence, 0.55)
    if (length(supports) >= 2L) {
        for (i in seq_len(length(supports) - 1L)) {
            if (i > boundary || !all(supports[i:(i + 1L)] < floor)) {
                next
            }
            boundary <- i - 1L
            break
        }
    }
    if (boundary < min_length) {
        return(0L)
    }
    called <- supports[seq_len(boundary)]
    if (mean(called >= confidence) < 0.8) {
        return(0L)
    }
    as.integer(boundary)
}
.detect_flank <- function(sequences, is_prefix, max_length = 60L, min_length = 14L,
    confidence = 0.75, min_sequences = 500L) {
    detected <- .positional_consensus(sequences, is_prefix, max_length, min_sequences)
    primer_length <- .flank_boundary(detected$supports, min_length, confidence)
    if (!primer_length) {
        return(list(sequence = NULL, length = 0L, confidence = 0))
    }
    bases <- detected$consensus[seq_len(primer_length)]
    primer <- paste0(if (is_prefix)
        bases else rev(bases), collapse = "")
    widths <- nchar(sequences)
    fragments <- if (is_prefix) {
        substr(sequences[widths >= primer_length], 1L, primer_length)
    } else {
        substr(sequences[widths >= primer_length], widths[widths >= primer_length] -
            primer_length + 1L, widths[widths >= primer_length])
    }
    if (length(fragments) < min_sequences) {
        return(list(sequence = NULL, length = 0L, confidence = 0))
    }
    list(sequence = primer, length = as.integer(primer_length), confidence = mean(vapply(fragments,
        .hamming_leq_one, logical(1), right = primer)))
}
.detect_primers <- function(sequences, confidence = 0.75, max_length = 60L, min_length = 14L,
    min_sequences = 500L) {
    if (length(sequences) < min_sequences) {
        return(NULL)
    }
    primer_5p <- .detect_flank(sequences, TRUE, max_length, min_length, confidence,
        min_sequences)
    primer_3p <- .detect_flank(sequences, FALSE, max_length, min_length, confidence,
        min_sequences)
    list(primer_5p = primer_5p, primer_3p = primer_3p, n_sequences_analyzed = length(sequences),
        mean_sequence_length = as.integer(round(mean(nchar(sequences)))), estimated_random_region_length = if (is.null(primer_5p$sequence) &&
            is.null(primer_3p$sequence)) NULL else as.integer(round(mean(nchar(sequences))) -
            primer_5p$length - primer_3p$length))
}
.normalise_pool <- function(sequences) {
    toupper(chartr("U", "T", sequences))
}
.adapter_reverse_complements <- function() {
    stats::setNames(vapply(.known_adapters, reverse_complement, character(1)),
        names(.known_adapters))
}
.matches_known_adapter_prefix <- function(primer) {
    if (is.null(primer) || !nzchar(primer)) {
        return(FALSE)
    }
    probe <- substr(primer, 1L, .adapter_probe_length)
    reverse_adapters <- .adapter_reverse_complements()
    any(vapply(.known_adapters, function(adapter) substr(adapter, 1L, .adapter_probe_length) ==
        probe, logical(1))) || any(vapply(reverse_adapters, function(adapter) {
        substr(adapter, 1L, .adapter_probe_length) == probe
    }, logical(1)))
}
.count_adapter_hits <- function(sequences) {
    reverse_adapters <- .adapter_reverse_complements()
    hits <- stats::setNames(vapply(names(.known_adapters), function(name) {
        forward <- substr(.known_adapters[[name]], 1L, .adapter_probe_length)
        reverse <- substr(reverse_adapters[[name]], 1L, .adapter_probe_length)
        sum(grepl(forward, sequences, fixed = TRUE) | grepl(reverse, sequences,
            fixed = TRUE))
    }, numeric(1)), names(.known_adapters))
    hits[order(names(hits), method = "radix")]
}
.top_variants <- function(sequences, length, is_prefix, k = 3L) {
    if (!length) {
        return(list())
    }
    widths <- nchar(sequences)
    fragments <- if (is_prefix) {
        substr(sequences[widths >= length], 1L, length)
    } else {
        substr(sequences[widths >= length], widths[widths >= length] - length +
            1L, widths[widths >= length])
    }
    if (!length(fragments)) {
        return(list())
    }
    values <- unique(fragments)
    counts <- tabulate(match(fragments, values), nbins = length(values))
    ordering <- order(-counts, seq_along(values), method = "radix")
    lapply(utils::head(ordering, k), function(index) list(values[[index]], as.numeric(counts[[index]])))
}
.position_consistency <- function(sequences, primer, is_prefix, tolerance = 3L) {
    if (is.null(primer)) {
        return(0)
    }
    primer_length <- nchar(primer)
    usable <- sequences[nchar(sequences) >= primer_length]
    if (!length(usable)) {
        return(0)
    }
    hits <- vapply(usable, function(sequence) {
        sequence_length <- nchar(sequence)
        for (offset in 0L:tolerance) {
            fragment <- if (is_prefix) {
                if (offset + primer_length > sequence_length)
                  next
                substr(sequence, offset + 1L, offset + primer_length)
            } else {
                if (primer_length + offset > sequence_length)
                  next
                substr(sequence, sequence_length - primer_length - offset + 1L,
                  sequence_length - offset)
            }
            if (.hamming_leq_one(fragment, primer)) {
                return(TRUE)
            }
        }
        FALSE
    }, logical(1))
    mean(hits)
}
.substring_match_rate <- function(sequences, primer) {
    if (is.null(primer)) {
        return(0)
    }
    primer_length <- nchar(primer)
    usable <- sequences[nchar(sequences) >= primer_length]
    if (!length(usable)) {
        return(0)
    }
    mean(vapply(usable, function(sequence) {
        sequence_length <- nchar(sequence)
        any(vapply(seq_len(sequence_length - primer_length + 1L), function(start) {
            .hamming_leq_one(substr(sequence, start, start + primer_length - 1L),
                primer)
        }, logical(1)))
    }, logical(1)))
}
.contains_hamming_leq_one <- function(longer, shorter) {
    if (nchar(shorter) > nchar(longer)) {
        return(FALSE)
    }
    any(vapply(seq_len(nchar(longer) - nchar(shorter) + 1L), function(start) {
        .hamming_leq_one(substr(longer, start, start + nchar(shorter) - 1L), shorter)
    }, logical(1)))
}
.primers_agree <- function(left, right) {
    if (is.null(left) || is.null(right)) {
        return(FALSE)
    }
    if (nchar(left) == nchar(right)) {
        return(.hamming_leq_one(left, right))
    }
    .contains_hamming_leq_one(if (nchar(left) > nchar(right))
        left else right, if (nchar(left) > nchar(right))
        right else left)
}
.persistence_score <- function(rates) {
    if (length(rates) < 2L || mean(rates) < 0.1) {
        return(NULL)
    }
    max(0, min(1, 1 - stats::sd(rates)/mean(rates)))
}
.combine_persistence <- function(persistence_5p, persistence_3p, primer_5p, primer_3p) {
    if (is.null(primer_5p) && is.null(primer_3p)) {
        return(NULL)
    }
    if (is.null(primer_5p)) {
        return(persistence_3p)
    }
    if (is.null(primer_3p)) {
        return(persistence_5p)
    }
    if (is.null(persistence_5p)) {
        return(persistence_3p)
    }
    if (is.null(persistence_3p)) {
        return(persistence_5p)
    }
    (persistence_5p + persistence_3p)/2
}
.n_length_stats <- function(sequences, primer_5p_length, primer_3p_length) {
    lengths <- pmax(0, nchar(sequences) - primer_5p_length - primer_3p_length)
    if (!length(lengths)) {
        return(list(mode = NULL, distribution = numeric(), confidence = 0))
    }
    values <- unique(lengths)
    counts <- tabulate(match(lengths, values), nbins = length(values))
    index <- which.max(counts)
    list(mode = as.numeric(values[[index]]), distribution = stats::setNames(as.numeric(counts),
        as.character(values)), confidence = counts[[index]]/sum(counts))
}
.detect_orientation <- function(sequences, primer_5p, primer_3p) {
    if (is.null(primer_5p) && is.null(primer_3p)) {
        return("FORWARD")
    }
    reverse_3p <- if (is.null(primer_3p))
        NULL else reverse_complement(primer_3p)
    forward <- 0L
    reverse <- 0L
    for (sequence in sequences) {
        if (!is.null(primer_5p) && nchar(sequence) >= nchar(primer_5p) && .hamming_leq_one(substr(sequence,
            1L, nchar(primer_5p)), primer_5p)) {
            forward <- forward + 1L
        } else if (!is.null(reverse_3p) && nchar(sequence) >= nchar(reverse_3p) &&
            .hamming_leq_one(substr(sequence, 1L, nchar(reverse_3p)), reverse_3p)) {
            reverse <- reverse + 1L
        }
    }
    total <- forward + reverse
    if (!total || reverse/total < 0.05) {
        return("FORWARD")
    }
    if (reverse/total > 0.95)
        "REVERSE" else "MIXED"
}
.composite_confidence <- function(signals, has_round_map) {
    weights <- if (has_round_map) {
        c(match_5p = 0.15, match_3p = 0.15, pos_5p = 0.15, pos_3p = 0.15, persistence = 0.25,
            n_length = 0.1, adapter_clean = 0.05)
    } else {
        c(match_5p = 0.225, match_3p = 0.225, pos_5p = 0.225, pos_3p = 0.225, persistence = 0,
            n_length = 0.05, adapter_clean = 0.05)
    }
    # Accumulate in schema order to match the reference implementation's
    # floating-point evaluation and reproducible JSON representation.
    total <- 0
    for (name in names(weights)) {
        value <- signals[[name]]
        if (!is.null(value)) {
            total <- total + weights[[name]] * value
        }
    }
    max(0, min(1, total))
}
.assign_status <- function(extraction_mode, confidence, has_round_map) {
    if (identical(extraction_mode, "UNABLE_TO_EXTRACT")) {
        return("UNABLE_TO_INFER")
    }
    status <- if (confidence >= 0.85) {
        "HIGH"
    } else if (confidence >= 0.6) {
        "MEDIUM"
    } else if (confidence >= 0.3) {
        "LOW"
    } else {
        "UNABLE_TO_INFER"
    }
    if (!has_round_map && identical(status, "HIGH"))
        "MEDIUM" else status
}
.classify_library <- function(match_rate_5p, match_rate_3p, n_length_confidence,
    has_paired_split, has_round_map, confidence) {
    if (match_rate_5p < 0.4 && match_rate_3p < 0.4) {
        return(list(extraction_mode = "UNABLE_TO_EXTRACT", full_insert_recovered = FALSE,
            required_action = "MANUAL_PRIMERS_REQUIRED", status = "UNABLE_TO_INFER",
            failure_reason = sprintf("Both primer match rates below 0.40 (5'=%.2f, 3'=%.2f)",
                match_rate_5p, match_rate_3p)))
    }
    if (has_paired_split) {
        mode <- "PAIRED_END_SPLIT_PRIMERS"
        return(list(extraction_mode = mode, full_insert_recovered = FALSE, required_action = "READ_MERGING_RECOMMENDED",
            status = .assign_status(mode, confidence, has_round_map), failure_reason = NULL))
    }
    if (match_rate_5p > 0.7 && match_rate_3p > 0.7) {
        mode <- "BOTH_PRIMERS_SINGLE_READ"
        return(list(extraction_mode = mode, full_insert_recovered = TRUE, required_action = "NONE",
            status = .assign_status(mode, confidence, has_round_map), failure_reason = NULL))
    }
    if (match_rate_5p > 0.7 || match_rate_3p > 0.7) {
        if (n_length_confidence > 0.8) {
            mode <- if (match_rate_5p > 0.7)
                "FIVE_PRIME_ONLY" else "THREE_PRIME_ONLY"
            return(list(extraction_mode = mode, full_insert_recovered = FALSE,
                required_action = "NONE", status = .assign_status(mode, confidence,
                  has_round_map), failure_reason = NULL))
        }
        side <- if (match_rate_5p > 0.7)
            "5'" else "3'"
        return(list(extraction_mode = "UNABLE_TO_EXTRACT", full_insert_recovered = FALSE,
            required_action = "MANUAL_PRIMERS_REQUIRED", status = "UNABLE_TO_INFER",
            failure_reason = sprintf("Only %s primer detected but N-length confidence (%.2f) <= threshold (0.80)",
                side, n_length_confidence)))
    }
    list(extraction_mode = "UNABLE_TO_EXTRACT", full_insert_recovered = FALSE,
        required_action = "MANUAL_PRIMERS_REQUIRED", status = "UNABLE_TO_INFER",
        failure_reason = sprintf("Ambiguous primer evidence: 5'=%.2f, 3'=%.2f, neither side passes the 0.70 'primer found' threshold",
            match_rate_5p, match_rate_3p))
}
.round_order <- function(round_names) {
    numeric_suffix <- suppressWarnings(as.numeric(sub(".*?([0-9]+)$", "\\1", round_names)))
    order(is.na(numeric_suffix), numeric_suffix, round_names, method = "radix")
}
.as_round_pools <- function(sequences_by_round, argument) {
    if (!is.list(sequences_by_round)) {
        stop(sprintf("`%s` must be a list of per-round sequence pools.", argument),
            call. = FALSE)
    }
    if (!length(sequences_by_round)) {
        return(list())
    }
    round_names <- names(sequences_by_round)
    if (is.null(round_names)) {
        round_names <- sprintf("round_%02d", seq_along(sequences_by_round) - 1L)
        names(sequences_by_round) <- round_names
    }
    if (any(!nzchar(round_names)) || anyDuplicated(round_names)) {
        stop(sprintf("`%s` must have unique, non-empty round names.", argument),
            call. = FALSE)
    }
    pools <- lapply(sequences_by_round, .as_sequence_character)
    pools[.round_order(names(pools))]
}
.with_sampling_seed <- function(seed, code) {
    withr::with_seed(as.integer(seed), code)
}
.subsample_pool <- function(sequences, max_reads) {
    if (is.null(max_reads) || max_reads <= 0L || length(sequences) <= max_reads) {
        return(list(sequences = sequences, fraction = 1))
    }
    list(sequences = sequences[sample.int(length(sequences), max_reads)], fraction = max_reads/length(sequences))
}
.paired_split_signal <- function(r1_sequences, r2_sequences) {
    r1_detection <- .detect_primers(r1_sequences)
    r2_detection <- .detect_primers(r2_sequences)
    if (is.null(r1_detection) || is.null(r2_detection)) {
        return(list(has_paired_split = FALSE, primer_3p = NULL))
    }
    r1_5p_strong <- r1_detection$primer_5p$confidence > 0.7
    r1_3p_strong <- r1_detection$primer_3p$confidence > 0.7
    r2_5p_strong <- r2_detection$primer_5p$confidence > 0.7
    r2_5p <- r2_detection$primer_5p$sequence
    if (!r1_5p_strong || !r2_5p_strong || is.null(r2_5p)) {
        return(list(has_paired_split = FALSE, primer_3p = NULL))
    }
    primer_3p <- reverse_complement(r2_5p)
    agrees <- .primers_agree(r1_detection$primer_3p$sequence, primer_3p)
    if (!r1_3p_strong || !agrees) {
        return(list(has_paired_split = TRUE, primer_3p = primer_3p))
    }
    list(has_paired_split = FALSE, primer_3p = NULL)
}
#' Detect SELEX library structure and constant regions
#'
#' Detects 5' and 3' constant regions in the earliest supplied SELEX round,
#' then verifies the call across rounds. The returned LibraryReport separates
#' inferred biology, input layout, and required downstream action. When primer
#' evidence is insufficient, it returns an explicit safe-failure report rather
#' than guessing.
#'
#' @param sequences_by_round A named list of character vectors or
#'   `Biostrings::XStringSet` objects. Each element is one sequencing round.
#' @param read_source One of `'R1'`, `'R2'`, `'R1_AND_R2'`, `'INTERLEAVED'`,
#'   or `'UNKNOWN'`.
#' @param paired_mate_streams Optional named list containing the corresponding
#'   R2 sequences. Required only to detect split-primer paired-end layouts.
#' @param sampling_seed Integer seed used when `max_reads_per_round` samples
#'   input reads.
#' @param max_reads_per_round Optional positive cap on reads used per round.
#'
#' @return A `selexprep_library_report`.
#' @export
#'
#' @examples
#' p5 <- 'GGTAATACGACTCACTATAGGG'
#' p3 <- 'CCATGCATGCATGCATGCAT'
#' selexprep_detect(list(round_00 = rep(paste0(p5, 'ACGTACGTACGTACGT', p3), 500)))
selexprep_detect <- function(sequences_by_round, read_source = c("R1", "R2", "R1_AND_R2",
    "INTERLEAVED", "UNKNOWN"), paired_mate_streams = NULL, sampling_seed = 42L,
    max_reads_per_round = NULL) {
    read_source <- match.arg(read_source)
    if (!is.numeric(sampling_seed) || length(sampling_seed) != 1L || is.na(sampling_seed)) {
        stop("`sampling_seed` must be one non-missing number.", call. = FALSE)
    }
    if (!is.null(max_reads_per_round) && (!is.numeric(max_reads_per_round) || length(max_reads_per_round) !=
        1L || is.na(max_reads_per_round))) {
        stop("`max_reads_per_round` must be NULL or one non-missing number.", call. = FALSE)
    }
    pools <- .as_round_pools(sequences_by_round, "sequences_by_round")
    if (!length(pools)) {
        return(.unable_library_report(read_source, sampling_seed, "No sequences provided"))
    }
    mate_pools <- if (is.null(paired_mate_streams))
        NULL else {
        .as_round_pools(paired_mate_streams, "paired_mate_streams")
    }
    if (!is.null(mate_pools) && !identical(names(mate_pools), names(pools))) {
        stop("`paired_mate_streams` must have the same round names as `sequences_by_round`.",
            call. = FALSE)
    }
    .with_sampling_seed(sampling_seed, {
        sampled_pools <- lapply(pools, function(pool) {
            .subsample_pool(.normalise_pool(pool), max_reads_per_round)
        })
        fractions <- vapply(sampled_pools, `[[`, numeric(1), "fraction")
        normalised <- lapply(sampled_pools, `[[`, "sequences")
        names(normalised) <- names(pools)
        earliest_name <- names(normalised)[[1L]]
        earliest_sequences <- normalised[[earliest_name]]
        adapter_hits <- .count_adapter_hits(earliest_sequences)
        detection <- .detect_primers(earliest_sequences)
        if (is.null(detection)) {
            failure_reason <- paste(
                "Earliest round has fewer than 500 sequences",
                "\u2014 below detection floor"
            )
            return(.unable_library_report(
                read_source,
                sampling_seed,
                failure_reason,
                adapter_hits
            ))
        }
        primer_5p <- detection$primer_5p$sequence
        primer_3p <- detection$primer_3p$sequence
        adapter_drop_5p <- .matches_known_adapter_prefix(primer_5p)
        adapter_drop_3p <- .matches_known_adapter_prefix(primer_3p)
        if (adapter_drop_5p)
            primer_5p <- NULL
        if (adapter_drop_3p)
            primer_3p <- NULL
        has_paired_split <- FALSE
        threep_sequences <- normalised
        threep_is_prefix <- FALSE
        primer_3p_lookup <- primer_3p
        if (!is.null(mate_pools)) {
            normalised_mates <- lapply(mate_pools, function(pool) {
                .subsample_pool(.normalise_pool(pool), max_reads_per_round)$sequences
            })
            split <- .paired_split_signal(earliest_sequences, normalised_mates[[earliest_name]])
            has_paired_split <- split$has_paired_split
            if (has_paired_split) {
                primer_3p <- split$primer_3p
                primer_3p_lookup <- reverse_complement(primer_3p)
                threep_sequences <- normalised_mates
                threep_is_prefix <- TRUE
            }
        }
        position_rates_5p <- vapply(normalised, .position_consistency, numeric(1),
            primer = primer_5p, is_prefix = TRUE)
        position_rates_3p <- vapply(threep_sequences, .position_consistency, numeric(1),
            primer = primer_3p_lookup, is_prefix = threep_is_prefix)
        persistence <- .combine_persistence(.persistence_score(position_rates_5p),
            .persistence_score(position_rates_3p), primer_5p, primer_3p)
        earliest_threep <- threep_sequences[[earliest_name]]
        match_rate_5p <- .substring_match_rate(earliest_sequences, primer_5p)
        match_rate_3p <- .substring_match_rate(earliest_threep, primer_3p_lookup)
        p5_length <- if (is.null(primer_5p))
            0L else nchar(primer_5p)
        p3_length <- if (is.null(primer_3p))
            0L else nchar(primer_3p)
        p3_lookup_length <- if (is.null(primer_3p_lookup))
            0L else nchar(primer_3p_lookup)
        variants_5p <- .top_variants(earliest_sequences, p5_length, TRUE)
        variants_3p <- .top_variants(earliest_threep, p3_lookup_length, threep_is_prefix)
        n_lengths <- if (has_paired_split) {
            list(mode = NULL, distribution = numeric(), confidence = 0)
        } else {
            .n_length_stats(earliest_sequences, p5_length, p3_length)
        }
        has_round_map <- length(normalised) >= 2L
        confidence <- .composite_confidence(list(match_5p = match_rate_5p, match_3p = match_rate_3p,
            pos_5p = position_rates_5p[[1L]], pos_3p = position_rates_3p[[1L]],
            persistence = persistence, n_length = n_lengths$confidence, adapter_clean = if (adapter_drop_5p ||
                adapter_drop_3p) 0 else 1), has_round_map)
        classification <- .classify_library(match_rate_5p, match_rate_3p, n_lengths$confidence,
            has_paired_split, has_round_map, confidence)
        .new_library_report(list(primer_5p = primer_5p, primer_3p = primer_3p,
            variants_5p = variants_5p, variants_3p = variants_3p, known_adapter_hits = adapter_hits,
            extraction_mode = classification$extraction_mode, full_insert_recovered = classification$full_insert_recovered,
            read_source = read_source, required_action = classification$required_action,
            orientation = .detect_orientation(earliest_sequences, primer_5p, primer_3p),
            n_length_mode = n_lengths$mode, n_length_distribution = n_lengths$distribution,
            n_length_confidence = n_lengths$confidence, match_rate_5p = match_rate_5p,
            match_rate_3p = match_rate_3p, position_consistency_5p = position_rates_5p[[1L]],
            position_consistency_3p = position_rates_3p[[1L]], read_fraction_used_for_inference = mean(fractions),
            sampling_seed = as.numeric(sampling_seed), confidence = confidence,
            status = classification$status, failure_reason = classification$failure_reason))
    })
}
.library_report_payload <- function(report) {
    if (!inherits(report, "selexprep_library_report")) {
        stop("`report` must be a selexprep_library_report.", call. = FALSE)
    }
    payload <- unclass(report)
    if (length(payload$n_length_distribution)) {
        ordering <- order(as.numeric(names(payload$n_length_distribution)), method = "radix")
        payload$n_length_distribution <- payload$n_length_distribution[ordering]
    }
    payload$known_adapter_hits <- payload$known_adapter_hits[order(names(payload$known_adapter_hits),
        method = "radix")]
    as_json_object <- function(values) {
        object <- lapply(unname(values), function(value) unname(value))
        names(object) <- names(values)
        object
    }
    payload$n_length_distribution <- as_json_object(payload$n_length_distribution)
    payload$known_adapter_hits <- as_json_object(payload$known_adapter_hits)
    payload
}
#' Write a LibraryReport as deterministic JSON
#'
#' @param report A `selexprep_library_report`.
#' @param path Output JSON path.
#'
#' @return The SHA-256 digest of the emitted UTF-8 JSON, invisibly.
#' @export
#' @examples
#' p5 <- 'GGTAATACGACTCACTATAGGG'
#' p3 <- 'CCATGCATGCATGCATGCAT'
#' report <- selexprep_detect(list(round_00 = rep(paste0(p5, 'ACGTACGTACGTACGT', p3), 500)))
#' path <- tempfile(fileext = '.json')
#' write_library_report(report, path)
write_library_report <- function(report, path) {
    payload <- .library_report_payload(report)
    json <- jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE, null = "null",
        digits = NA)
    # jsonlite represents an empty named list as [] even though these stable
    # schema fields are maps. Preserve the Python-compatible empty object.
    json <- sub("\"known_adapter_hits\": []", "\"known_adapter_hits\": {}", json,
        fixed = TRUE)
    json <- sub("\"n_length_distribution\": []", "\"n_length_distribution\": {}",
        json, fixed = TRUE)
    # Python's schema writes float-valued signals as `1.0` rather than `1`.
    # Keep that stable spelling and its 16-significant-digit representation
    # so JSON produced by both implementations is byte-identical for
    # reference reports.
    float_fields <- c("n_length_confidence", "match_rate_5p", "match_rate_3p",
        "position_consistency_5p", "position_consistency_3p", "read_fraction_used_for_inference",
        "confidence")
    format_python_float <- function(value) {
        formatted <- sprintf("%.16g", value)
        if (!grepl("[.eE]", formatted)) {
            formatted <- paste0(formatted, ".0")
        }
        formatted
    }
    for (field in float_fields) {
        pattern <- sprintf("(\\\"%s\\\":\\s*)-?[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?",
            field)
        json <- sub(pattern, paste0("\\1", format_python_float(payload[[field]])),
            json, perl = TRUE)
    }
    text <- paste0(json, "\n")
    .write_utf8_lf(text, path)
    invisible(digest::digest(text, algo = "sha256", serialize = FALSE))
}
#' Read a LibraryReport JSON file
#'
#' @param path Path to a JSON file written by `write_library_report()` or by
#'   the compatible Python implementation.
#'
#' @return A validated `selexprep_library_report`.
#' @export
#' @examples
#' p5 <- 'GGTAATACGACTCACTATAGGG'
#' p3 <- 'CCATGCATGCATGCATGCAT'
#' report <- selexprep_detect(list(round_00 = rep(paste0(p5, 'ACGTACGTACGTACGT', p3), 500)))
#' path <- tempfile(fileext = '.json')
#' write_library_report(report, path)
#' read_library_report(path)
read_library_report <- function(path) {
    payload <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    required <- .library_report_fields
    if (!identical(names(payload), required)) {
        stop("The JSON file does not match the stable LibraryReport schema.", call. = FALSE)
    }
    as_named_numeric <- function(values) {
        flattened <- unlist(values, use.names = TRUE)
        stats::setNames(as.numeric(flattened), names(flattened))
    }
    # `[<-` preserves a zero-length field, whereas `$<-` would remove it.
    payload["known_adapter_hits"] <- list(as_named_numeric(payload$known_adapter_hits))
    payload["n_length_distribution"] <- list(as_named_numeric(payload$n_length_distribution))
    .new_library_report(payload)
}
