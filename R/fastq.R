#' Read sequences from a FASTQ file
#'
#' Reads a plain or gzip-compressed FASTQ file into a `Biostrings::BStringSet`.
#' `BStringSet` deliberately preserves both DNA and RNA alphabets, avoiding the
#' loss of uracil that would occur when every input is coerced to DNA.
#'
#' @param path Path to a FASTQ or FASTQ.gz file.
#' @param max_reads Optional positive maximum number of records to retain.
#'
#' @return A `Biostrings::BStringSet` of read sequences.
#' @export
#' @examples
#' path <- tempfile(fileext = '.fastq')
#' writeLines(c('@read_1', 'ACGU', '+', 'IIII'), path)
#' read_selexprep_fastq(path)
read_selexprep_fastq <- function(path, max_reads = NULL) {
    valid_path <- is.character(path) && length(path) == 1L && !is.na(path) && file.exists(path)
    if (!valid_path) {
        stop("`path` must name an existing FASTQ file.", call. = FALSE)
    }
    if (!is.null(max_reads) && (!is.numeric(max_reads) || length(max_reads) !=
        1L || is.na(max_reads) || max_reads < 1L || max_reads != floor(max_reads))) {
        stop("`max_reads` must be NULL or one positive integer.", call. = FALSE)
    }
    nrec <- if (is.null(max_reads))
        -1L else as.integer(max_reads)
    Biostrings::readBStringSet(path, format = "fastq", nrec = nrec)
}
