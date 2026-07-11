# Convert supported sequence containers to un-named character vectors.
.as_sequence_character <- function(sequences) {
    if (is.character(sequences)) {
        values <- unname(sequences)
    } else if (inherits(sequences, "XStringSet")) {
        values <- unname(as.character(sequences))
    } else {
        stop("`sequences` must be a character vector or a Biostrings XStringSet.",
            call. = FALSE)
    }
    if (anyNA(values)) {
        stop("`sequences` must not contain missing values.", call. = FALSE)
    }
    values
}
#' Count unique SELEX sequences
#'
#' Counts the observed sequences in one primer-stripped SELEX round. The
#' result is a `S4Vectors::DataFrame` compatible with the stable Python output
#' schema: `sequence`, `reads`, `rank`, and `rpm`. The sequence column is a
#' `Biostrings::BStringSet` so that both DNA and RNA alphabets can be carried
#' without lossy coercion.
#'
#' Ties in abundance are ordered lexicographically by sequence. The original
#' Python implementation does not define a tie-breaker; making it explicit in
#' R provides reproducible ranks across platforms.
#'
#' @param sequences A character vector or a `Biostrings::XStringSet` containing
#'   the primer-stripped sequences from one round.
#'
#' @return A `S4Vectors::DataFrame` with `sequence`, `reads`, `rank`, and `rpm`
#'   columns, ordered by decreasing read count.
#' @export
#'
#' @examples
#' selexprep_count(c('ACGT', 'ACGT', 'GGGG'))
selexprep_count <- function(sequences) {
    values <- .as_sequence_character(sequences)
    if (!length(values)) {
        return(S4Vectors::DataFrame(sequence = Biostrings::BStringSet(), reads = numeric(),
            rank = integer(), rpm = numeric()))
    }
    unique_values <- unique(values)
    reads <- tabulate(match(values, unique_values), nbins = length(unique_values))
    ordering <- order(-reads, unique_values, method = "radix")
    unique_values <- unique_values[ordering]
    reads <- as.numeric(reads[ordering])
    n_reads <- sum(reads)
    S4Vectors::DataFrame(sequence = Biostrings::BStringSet(unique_values), reads = reads,
        rank = seq_along(reads), rpm = reads/n_reads * 1e+06)
}
# Assemble one count table per round into the standard multi-round container.
.as_selexprep_experiment <- function(counts_by_round, accession = NULL) {
    if (!is.list(counts_by_round) || !length(counts_by_round)) {
        stop("`counts_by_round` must be a non-empty named list of count tables.",
            call. = FALSE)
    }
    round_names <- names(counts_by_round)
    if (is.null(round_names) || any(!nzchar(round_names)) || anyDuplicated(round_names)) {
        stop("`counts_by_round` must have unique, non-empty round names.", call. = FALSE)
    }
    sequences_by_round <- lapply(counts_by_round, function(counts) {
        if (!inherits(counts, "DataFrame") || !all(c("sequence", "reads") %in%
            colnames(counts))) {
            stop("Each element of `counts_by_round` must be a DataFrame from `selexprep_count()`.",
                call. = FALSE)
        }
        as.character(counts$sequence)
    })
    all_sequences <- sort(unique(unlist(sequences_by_round, use.names = FALSE)))
    row_index <- unlist(lapply(sequences_by_round, match, table = all_sequences),
        use.names = FALSE)
    col_index <- rep(seq_along(counts_by_round), lengths(sequences_by_round))
    count_values <- unlist(lapply(counts_by_round, function(counts) as.numeric(counts$reads)),
        use.names = FALSE)
    counts <- Matrix::sparseMatrix(i = row_index, j = col_index, x = count_values,
        dims = c(length(all_sequences), length(counts_by_round)), dimnames = list(all_sequences,
            round_names))
    SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts),
        rowData = S4Vectors::DataFrame(sequence = Biostrings::BStringSet(all_sequences),
            row.names = all_sequences), colData = S4Vectors::DataFrame(round = round_names,
            row.names = round_names), metadata = list(accession = accession, schema_version = "0.1.0"))
}
