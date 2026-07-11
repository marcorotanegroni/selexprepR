.fastq_mate_from_path <- function(path) {
    stem <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(path), ignore.case = TRUE)
    marked <- regexec("(?:^|[._-])R([12])(?:[._-]|$)", stem, ignore.case = TRUE,
        perl = TRUE)
    value <- regmatches(stem, marked)[[1L]]
    if (length(value) >= 2L) {
        return(paste0("R", value[[2L]]))
    }
    trailing <- regexec("(?:^|[._-])([12])$", stem, perl = TRUE)
    value <- regmatches(stem, trailing)[[1L]]
    if (length(value) >= 2L)
        paste0("R", value[[2L]]) else NA_character_
}
.normalise_fastq_mates <- function(mate, paths) {
    if (is.null(mate)) {
        mate <- vapply(paths, .fastq_mate_from_path, character(1))
    } else {
        if (length(mate) == 1L) {
            mate <- rep(mate, length(paths))
        }
        if (length(mate) != length(paths) || anyNA(mate)) {
            stop("`mate` must contain one value per FASTQ file.", call. = FALSE)
        }
        mate <- toupper(as.character(mate))
        mate[mate == "1"] <- "R1"
        mate[mate == "2"] <- "R2"
        if (any(!mate %in% c("R1", "R2"))) {
            stop("`mate` values must be R1, R2, 1, or 2.", call. = FALSE)
        }
    }
    if (anyNA(mate)) {
        if (any(mate == "R2", na.rm = TRUE)) {
            stop("Could not infer every mate in a paired-end input; supply `mate` explicitly.",
                call. = FALSE)
        }
        mate[is.na(mate)] <- "R1"
    }
    mate
}
.validate_fastq_rounds <- function(round, n_files) {
    if (length(round) == 1L) {
        round <- rep(round, n_files)
    }
    valid <- is.numeric(round) && length(round) == n_files && !anyNA(round) &&
        all(is.finite(round)) && all(round >= 0) && all(round == floor(round))
    if (!valid) {
        stop("`round` must contain one non-negative integer per FASTQ file.", call. = FALSE)
    }
    as.integer(round)
}
.rounds_from_path_names <- function(paths) {
    labels <- names(paths)
    if (is.null(labels) || any(!nzchar(labels))) {
        return(NULL)
    }
    if (!all(grepl("[0-9]+$", labels))) {
        return(NULL)
    }
    suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", labels)))
}
.fetch_result_input_table <- function(x) {
    if (isTRUE(x$dry_run)) {
        stop("A dry-run fetch result has no FASTQ files; download it first.", call. = FALSE)
    }
    plan <- x$plan
    paths <- x$downloaded_files
    valid_plan <- inherits(plan, "DataFrame") || is.data.frame(plan)
    required <- c("round_number", "file_name")
    if (!valid_plan || !all(required %in% colnames(plan))) {
        stop("The fetch result contains an invalid download plan.", call. = FALSE)
    }
    if (!is.character(paths) || length(paths) != nrow(plan)) {
        stop("The fetch result does not map every plan row to a downloaded file.",
            call. = FALSE)
    }
    if ("round_unambiguous" %in% colnames(plan)) {
        unambiguous <- as.logical(plan$round_unambiguous)
        if (anyNA(unambiguous) || any(!unambiguous)) {
            stop("The fetch result contains ambiguous round assignments.", call. = FALSE)
        }
    }
    run_accession <- if ("run_accession" %in% colnames(plan)) {
        as.character(plan$run_accession)
    } else {
        rep(NA_character_, nrow(plan))
    }
    data.frame(path = paths, round_number = plan$round_number, mate = vapply(as.character(plan$file_name),
        .fastq_mate_from_path, character(1)), run_accession = run_accession, stringsAsFactors = FALSE)
}
.local_input_table <- function(x, round, mate) {
    if (is.character(x)) {
        paths <- x
        if (is.null(round)) {
            round <- .rounds_from_path_names(paths)
        }
        if (is.null(round)) {
            stop("Supply `round`, or name each path with a label ending in its round number.",
                call. = FALSE)
        }
        return(data.frame(path = unname(paths), round_number = round, mate = if (is.null(mate)) NA_character_ else mate,
            run_accession = NA_character_, stringsAsFactors = FALSE))
    }
    valid_table <- inherits(x, "DataFrame") || is.data.frame(x)
    if (!valid_table) {
        stop("`x` must be FASTQ paths, a file table, or a completed selexprep_fetch_result.",
            call. = FALSE)
    }
    required <- c("path", "round_number")
    if (!all(required %in% colnames(x))) {
        stop("An input file table must contain `path` and `round_number`.", call. = FALSE)
    }
    n_files <- nrow(x)
    data.frame(path = as.character(x$path), round_number = if (is.null(round))
        x$round_number else round, mate = if (!is.null(mate)) {
        mate
    } else if ("mate" %in% colnames(x)) {
        as.character(x$mate)
    } else {
        rep(NA_character_, n_files)
    }, run_accession = if ("run_accession" %in% colnames(x)) {
        as.character(x$run_accession)
    } else {
        rep(NA_character_, n_files)
    }, stringsAsFactors = FALSE)
}
.combine_fastq_sets <- function(sets) {
    if (length(sets) == 1L) {
        return(sets[[1L]])
    }
    do.call(c, unname(sets))
}
#' Read FASTQ files into round-aware pipeline inputs
#'
#' Converts local FASTQ paths or a completed ENA fetch result into the R-native
#' inputs used by [run_selexprep()]. Files belonging to the same round and mate
#' are concatenated in input-table order. Plain and gzip-compressed FASTQ files
#' are supported through [read_selexprep_fastq()].
#'
#' Mate labels are inferred conservatively from conventional `_R1`, `_R2`,
#' `_1`, and `_2` filename suffixes. Files without a mate suffix are treated as
#' single-end only when no R2 file is present. Paired inputs must have both
#' mates, with equal loaded read counts, in every round.
#'
#' @param x A character vector of FASTQ paths, a data frame-like object with
#'   `path` and `round_number` columns, or a non-dry-run
#'   `selexprep_fetch_result` returned by [fetch_selexprep_reads()].
#' @param round Optional vector of non-negative round numbers. It is required
#'   for character paths unless their vector names end in round numbers.
#' @param mate Optional vector containing `R1`/`R2` or `1`/`2`. When omitted,
#'   mate labels are inferred from filenames.
#' @param max_reads_per_file Optional positive maximum number of records read
#'   from each FASTQ file.
#'
#' @return A `selexprep_read_inputs` list with `sequences_by_round`, optional
#'   `paired_mate_streams`, `read_source`, and a file-level `DataFrame` named
#'   `files`. The first three elements can be passed directly to
#'   [run_selexprep()], or use [run_selexprep_files()].
#' @export
#'
#' @examples
#' path <- tempfile(fileext = '.fastq')
#' writeLines(c('@read_1', 'ACGT', '+', 'IIII'), path)
#' inputs <- read_selexprep_inputs(path, round = 0)
#' inputs$sequences_by_round
read_selexprep_inputs <- function(x, round = NULL, mate = NULL, max_reads_per_file = NULL) {
    table <- if (inherits(x, "selexprep_fetch_result")) {
        .fetch_result_input_table(x)
    } else {
        .local_input_table(x, round = round, mate = mate)
    }
    n_files <- nrow(table)
    if (!n_files) {
        stop("At least one FASTQ file is required.", call. = FALSE)
    }
    table$round_number <- .validate_fastq_rounds(table$round_number, n_files)
    valid_paths <- is.character(table$path) && length(table$path) == n_files &&
        !anyNA(table$path) && all(nzchar(table$path)) && all(file.exists(table$path))
    if (!valid_paths) {
        stop("Every `path` must name an existing FASTQ file.", call. = FALSE)
    }
    table$path <- normalizePath(table$path, winslash = "/", mustWork = TRUE)
    if (anyDuplicated(table$path)) {
        stop("Each FASTQ path may appear only once.", call. = FALSE)
    }
    supplied_mate <- if (all(is.na(table$mate)))
        NULL else table$mate
    table$mate <- .normalise_fastq_mates(supplied_mate, table$path)
    paired <- any(table$mate == "R2")
    rounds <- sort(unique(table$round_number), method = "radix")
    if (paired) {
        complete <- vapply(rounds, function(value) {
            labels <- unique(table$mate[table$round_number == value])
            identical(sort(labels), c("R1", "R2"))
        }, logical(1))
        if (!all(complete)) {
            stop("Paired-end inputs require both R1 and R2 in every round.", call. = FALSE)
        }
    }
    sets <- lapply(table$path, read_selexprep_fastq, max_reads = max_reads_per_file)
    table$reads_loaded <- vapply(sets, length, integer(1))
    table$file_size <- unname(file.info(table$path)$size)
    table$md5 <- vapply(table$path, digest::digest, character(1), algo = "md5",
        file = TRUE)
    table$sha256 <- vapply(table$path, digest::digest, character(1), algo = "sha256",
        file = TRUE)
    round_names <- sprintf("round_%02d", rounds)
    combine_mate <- function(label) {
        pools <- lapply(rounds, function(value) {
            keep <- table$round_number == value & table$mate == label
            .combine_fastq_sets(sets[keep])
        })
        stats::setNames(pools, round_names)
    }
    sequences_by_round <- combine_mate("R1")
    paired_mate_streams <- if (paired)
        combine_mate("R2") else NULL
    if (paired) {
        r1_reads <- vapply(sequences_by_round, length, integer(1))
        r2_reads <- vapply(paired_mate_streams, length, integer(1))
        unequal <- r1_reads != r2_reads
        if (any(unequal)) {
            stop(sprintf("Paired streams have different read counts in %s.",
                toString(names(r1_reads)[unequal])), call. = FALSE)
        }
    }
    file_order <- order(table$round_number, match(table$mate, c("R1", "R2")), seq_len(n_files),
        method = "radix")
    files <- S4Vectors::DataFrame(table[file_order, , drop = FALSE])
    structure(list(sequences_by_round = sequences_by_round, paired_mate_streams = paired_mate_streams,
        read_source = if (paired) "R1_AND_R2" else "R1", files = files), class = c("selexprep_read_inputs",
        "list"))
}
#' Run selexprepR directly from FASTQ files
#'
#' Convenience bridge combining [read_selexprep_inputs()] and
#' [run_selexprep()]. File-level provenance is attached to
#' `metadata(result)$input_files`.
#'
#' @inheritParams read_selexprep_inputs
#' @param ... Named arguments forwarded to [run_selexprep()]. The input
#'   sequence, paired-mate, and read-source arguments are supplied by this
#'   function and cannot be overridden.
#'
#' @return A `SummarizedExperiment` as returned by [run_selexprep()].
#' @export
#'
#' @examples
#' \donttest{
#' experiment <- run_selexprep_files(
#'     c(round_00 = 'round_00.fastq.gz',
#'       round_01 = 'round_01.fastq.gz')
#' )
#' }
run_selexprep_files <- function(x, round = NULL, mate = NULL, max_reads_per_file = NULL,
    ...) {
    inputs <- read_selexprep_inputs(x, round = round, mate = mate, max_reads_per_file = max_reads_per_file)
    extras <- list(...)
    if (length(extras) && (is.null(names(extras)) || any(!nzchar(names(extras))))) {
        stop("Arguments in `...` must be named.", call. = FALSE)
    }
    reserved <- c("sequences_by_round", "paired_mate_streams", "read_source")
    if (length(intersect(names(extras), reserved))) {
        stop("Sequence layout arguments are derived from the FASTQ inputs.", call. = FALSE)
    }
    arguments <- c(list(sequences_by_round = inputs$sequences_by_round, read_source = inputs$read_source,
        paired_mate_streams = inputs$paired_mate_streams), extras)
    experiment <- do.call(run_selexprep, arguments)
    metadata <- S4Vectors::metadata(experiment)
    metadata$input_files <- inputs$files
    manifest_fields <- unclass(metadata$manifest)
    file_keys <- sprintf("round_%02d/%s/%03d_%s", inputs$files$round_number, inputs$files$mate,
        seq_len(nrow(inputs$files)), basename(inputs$files$path))
    manifest_fields$input_sha256 <- stats::setNames(as.character(inputs$files$sha256),
        file_keys)
    metadata$manifest <- .new_selexprep_manifest(manifest_fields)
    S4Vectors::metadata(experiment) <- metadata
    experiment
}
