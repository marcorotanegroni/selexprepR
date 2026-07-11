.validate_demultiplex_barcodes <- function(barcodes, max_mismatches) {
    valid_rounds <- is.numeric(barcodes) && length(barcodes) > 0L &&
        !is.null(names(barcodes)) && all(nzchar(names(barcodes))) &&
        !anyNA(barcodes) && all(is.finite(barcodes)) &&
        all(barcodes >= 0) && all(barcodes == floor(barcodes))
    if (!valid_rounds || anyDuplicated(names(barcodes))) {
        stop(
            "`barcodes` must be a named non-negative integer vector.",
            call. = FALSE
        )
    }
    valid_mismatches <- is.numeric(max_mismatches) &&
        length(max_mismatches) == 1L && !is.na(max_mismatches) &&
        max_mismatches >= 0 && max_mismatches == floor(max_mismatches)
    if (!valid_mismatches) {
        stop(
            "`max_mismatches` must be one non-negative integer.",
            call. = FALSE
        )
    }
    barcode_sequences <- toupper(chartr("U", "T", names(barcodes)))
    widths <- nchar(barcode_sequences)
    if (length(unique(widths)) != 1L || any(widths == 0L) ||
        any(grepl("[^ACGT]", barcode_sequences))) {
        stop(
            "Barcode names must be equal-length DNA or RNA sequences.",
            call. = FALSE
        )
    }
    minimum_distance <- 2L * as.integer(max_mismatches) + 1L
    if (length(barcode_sequences) > 1L) {
        for (left in seq_len(length(barcode_sequences) - 1L)) {
            for (right in seq.int(left + 1L, length(barcode_sequences))) {
                distance <- sum(
                    strsplit(barcode_sequences[[left]], "")[[1L]] !=
                        strsplit(barcode_sequences[[right]], "")[[1L]]
                )
                if (distance < minimum_distance) {
                    stop(
                        sprintf(
                            paste(
                                "Barcodes %s and %s have distance %d;",
                                "at least %d is required."
                            ),
                            names(barcodes)[[left]],
                            names(barcodes)[[right]],
                            distance,
                            minimum_distance
                        ),
                        call. = FALSE
                    )
                }
            }
        }
    }
    stats::setNames(as.integer(barcodes), barcode_sequences)
}

.demultiplex_assignments <- function(sequences, barcodes, max_mismatches) {
    barcode_width <- nchar(names(barcodes)[[1L]])
    vapply(sequences, function(sequence) {
        sequence <- toupper(chartr("U", "T", sequence))
        if (nchar(sequence) < barcode_width) {
            return(NA_integer_)
        }
        prefix <- strsplit(substr(sequence, 1L, barcode_width), "")[[1L]]
        distances <- vapply(names(barcodes), function(barcode) {
            sum(prefix != strsplit(barcode, "")[[1L]])
        }, integer(1))
        matched <- which(distances <= max_mismatches)
        if (length(matched) == 1L) barcodes[[matched]] else NA_integer_
    }, integer(1))
}

#' Demultiplex SELEX reads with user-supplied barcodes
#'
#' Assigns reads to selection rounds from an equal-length barcode at the R1
#' 5-prime end. Barcode sequences are supplied explicitly; they are never
#' inferred. R2 mates follow their R1 assignment and remain untrimmed.
#'
#' Barcodes must have pairwise Hamming distance of at least
#' `2 * max_mismatches + 1`, preventing ambiguous assignments. RNA `U` and DNA
#' `T` are treated as equivalent during matching.
#'
#' @param sequences R1 sequences as a character vector or
#'   `Biostrings::XStringSet`.
#' @param barcodes Named integer vector mapping barcode sequences to
#'   non-negative round numbers, for example `c(AAAAA = 0L, TTTTT = 1L)`.
#' @param paired_mates Optional R2 sequences in the same order as `sequences`.
#' @param max_mismatches Maximum Hamming distance allowed from a barcode.
#' @param trim_barcode Whether to remove the matched barcode from assigned R1
#'   sequences.
#'
#' @return A `selexprep_demultiplex` object. Its first four fields follow the
#'   `selexprep_read_inputs` contract and can be supplied directly to
#'   [run_selexprep()]. Additional fields contain unassigned reads, a summary,
#'   and the validated barcode map.
#' @export
#'
#' @examples
#' reads <- c("AAAAACCCCC", "AAAATGGGGG", "TTTTTACGTA", "GGGGGAAAAA")
#' demultiplexed <- selexprep_demultiplex(
#'     reads,
#'     c(AAAAA = 0L, TTTTT = 1L)
#' )
#' demultiplexed$summary
selexprep_demultiplex <- function(sequences, barcodes, paired_mates = NULL,
    max_mismatches = 1L, trim_barcode = TRUE) {
    sequences <- .as_sequence_character(sequences)
    if (!is.null(paired_mates)) {
        paired_mates <- .as_sequence_character(paired_mates)
        if (length(sequences) != length(paired_mates)) {
            stop(
                "`paired_mates` must contain one mate per R1 sequence.",
                call. = FALSE
            )
        }
    }
    if (!is.logical(trim_barcode) || length(trim_barcode) != 1L ||
        is.na(trim_barcode)) {
        stop("`trim_barcode` must be one logical value.", call. = FALSE)
    }
    barcodes <- .validate_demultiplex_barcodes(
        barcodes,
        max_mismatches
    )
    assignments <- .demultiplex_assignments(
        sequences,
        barcodes,
        as.integer(max_mismatches)
    )
    rounds <- sort(unique(unname(barcodes)), method = "radix")
    round_names <- sprintf("round_%02d", rounds)
    barcode_width <- nchar(names(barcodes)[[1L]])
    r1 <- lapply(rounds, function(round) {
        selected <- sequences[!is.na(assignments) & assignments == round]
        if (trim_barcode) {
            selected <- substring(selected, barcode_width + 1L)
        }
        Biostrings::BStringSet(unname(selected))
    })
    names(r1) <- round_names
    r2 <- if (is.null(paired_mates)) {
        NULL
    } else {
        pools <- lapply(rounds, function(round) {
            keep <- !is.na(assignments) & assignments == round
            Biostrings::BStringSet(unname(paired_mates[keep]))
        })
        stats::setNames(pools, round_names)
    }
    assigned <- vapply(rounds, function(round) {
        sum(assignments == round, na.rm = TRUE)
    }, integer(1))
    total <- length(sequences)
    summary <- S4Vectors::DataFrame(
        round = round_names,
        round_number = rounds,
        assigned_reads = assigned,
        assigned_fraction = if (total) {
            assigned / total
        } else {
            numeric(length(rounds))
        }
    )
    unassigned <- is.na(assignments)
    unassigned_reads <- list(
        r1 = Biostrings::BStringSet(unname(sequences[unassigned])),
        r2 = if (is.null(paired_mates)) {
            NULL
        } else {
            Biostrings::BStringSet(unname(paired_mates[unassigned]))
        }
    )
    structure(list(
        sequences_by_round = r1,
        paired_mate_streams = r2,
        read_source = if (is.null(r2)) "R1" else "R1_AND_R2",
        files = S4Vectors::DataFrame(),
        unassigned = unassigned_reads,
        summary = summary,
        barcodes = barcodes,
        max_mismatches = as.integer(max_mismatches),
        trim_barcode = trim_barcode
    ), class = c(
        "selexprep_demultiplex",
        "selexprep_read_inputs",
        "list"
    ))
}
