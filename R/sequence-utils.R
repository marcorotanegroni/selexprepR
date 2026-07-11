#' Reverse-complement nucleotide sequences
#'
#' This function reproduces the sequence-counting behavior of the Python
#' implementation. It accepts DNA or RNA bases, keeps input case, maps U/u to
#' A/a, and preserves ambiguous bases such as N/n. It is vectorized over
#' `sequence`.
#'
#' @param sequence A character vector of nucleotide sequences without missing
#'   values.
#'
#' @return A character vector of reverse-complemented sequences.
#' @export
#'
#' @examples
#' reverse_complement(c('ACGT', 'ACGU', 'ANNT'))
reverse_complement <- function(sequence) {
    if (!is.character(sequence)) {
        stop("`sequence` must be a character vector.", call. = FALSE)
    }
    if (anyNA(sequence)) {
        stop("`sequence` must not contain missing values.", call. = FALSE)
    }
    complements <- chartr("ACGTUNacgtun", "TGCAANtgcaan", sequence)
    vapply(strsplit(complements, "", fixed = TRUE), function(bases) paste0(rev(bases),
        collapse = ""), character(1))
}
