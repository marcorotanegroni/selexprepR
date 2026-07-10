.as_selexprep_count_table <- function(counts) {
    if (!inherits(counts, "DataFrame") && !is.data.frame(counts)) {
        stop("`counts` must be a DataFrame or data.frame.", call. = FALSE)
    }
    expected <- c("sequence", "reads", "rank", "rpm")
    if (!identical(colnames(counts), expected)) {
        stop("Counts must have the stable columns: sequence, reads, rank, rpm.", call. = FALSE)
    }
    sequence <- as.character(counts$sequence)
    reads <- as.numeric(counts$reads)
    rank <- as.numeric(counts$rank)
    rpm <- as.numeric(counts$rpm)
    if (anyNA(sequence) || anyNA(reads) || anyNA(rank) || anyNA(rpm) ||
        any(!is.finite(reads)) || any(!is.finite(rank)) || any(!is.finite(rpm)) ||
        any(reads < 0) || any(rank < 1) || any(rank != floor(rank))) {
        stop("Counts contain invalid sequence, read, rank, or RPM values.", call. = FALSE)
    }
    if (!identical(rank, seq_along(rank) * 1)) {
        stop("Count ranks must be consecutive integers starting at one.", call. = FALSE)
    }
    total_reads <- sum(reads)
    expected_rpm <- if (total_reads) reads / total_reads * 1e6 else numeric()
    if (length(rpm) && any(abs(rpm - expected_rpm) > 1e-7 * pmax(1, expected_rpm))) {
        stop("Count RPM values do not agree with the read totals.", call. = FALSE)
    }
    S4Vectors::DataFrame(
        sequence = Biostrings::BStringSet(sequence),
        reads = reads,
        rank = as.integer(rank),
        rpm = rpm
    )
}

.count_path_format <- function(path) {
    name <- tolower(basename(path))
    if (grepl("\\.tsv$", name)) return("tsv")
    if (grepl("\\.csv$", name)) return("csv")
    if (grepl("\\.rds$", name)) return("rds")
    if (grepl("\\.parquet$", name)) return("parquet")
    stop("Unsupported count-file extension; use .tsv, .csv, .rds, or .parquet.", call. = FALSE)
}

.require_arrow <- function() {
    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop(
            "Parquet support requires the optional `arrow` package.",
            call. = FALSE
        )
    }
}

#' Read a selexprep count table
#'
#' Reads the stable four-column count contract (`sequence`, `reads`, `rank`,
#' `rpm`) from TSV, CSV, RDS, or Parquet. Parquet support is optional and
#' requires `arrow`; TSV is the portable default for R-native workflows.
#'
#' @param path Path to a `.tsv`, `.csv`, `.rds`, or `.parquet` count file.
#'
#' @return A validated `S4Vectors::DataFrame` with a
#'   `Biostrings::BStringSet` sequence column.
#' @export
#' @examples
#' path <- tempfile(fileext = ".tsv")
#' write_selexprep_counts(selexprep_count(c("ACGT", "ACGT", "GGGG")), path)
#' read_selexprep_counts(path)
read_selexprep_counts <- function(path) {
    format <- .count_path_format(path)
    counts <- switch(format,
        tsv = utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE),
        csv = utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
        rds = readRDS(path),
        parquet = {
            .require_arrow()
            getExportedValue("arrow", "read_parquet")(path, as_data_frame = TRUE)
        }
    )
    .as_selexprep_count_table(counts)
}

#' Write a selexprep count table
#'
#' Writes the stable four-column count contract. Use TSV for a dependency-free,
#' portable exchange format; Parquet is available when the optional `arrow`
#' package is installed.
#'
#' @param counts A count table returned by `selexprep_count()`.
#' @param path Output path ending in `.tsv`, `.csv`, `.rds`, or `.parquet`.
#'
#' @return `path`, invisibly.
#' @export
#' @examples
#' counts <- selexprep_count(c("ACGT", "ACGT", "GGGG"))
#' path <- tempfile(fileext = ".tsv")
#' write_selexprep_counts(counts, path)
write_selexprep_counts <- function(counts, path) {
    counts <- .as_selexprep_count_table(counts)
    format <- .count_path_format(path)
    payload <- data.frame(
        sequence = as.character(counts$sequence),
        reads = as.numeric(counts$reads),
        rank = as.integer(counts$rank),
        rpm = as.numeric(counts$rpm),
        check.names = FALSE
    )
    switch(format,
        tsv = utils::write.table(payload, path, sep = "\t", quote = FALSE, row.names = FALSE),
        csv = utils::write.csv(payload, path, quote = FALSE, row.names = FALSE),
        rds = saveRDS(counts, path),
        parquet = {
            .require_arrow()
            getExportedValue("arrow", "write_parquet")(payload, path)
        }
    )
    invisible(path)
}
