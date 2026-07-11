.ena_filereport_url <- "https://www.ebi.ac.uk/ena/portal/api/filereport"
.ena_filereport_fields <- c("run_accession", "study_accession", "study_title",
    "library_strategy", "library_source", "library_name", "experiment_title", "sample_title",
    "sample_accession", "read_count", "base_count", "fastq_md5", "fastq_bytes",
    "fastq_ftp")
.as_ena_scalar <- function(value) {
    if (!length(value) || is.na(value))
        "" else as.character(value)
}
.ena_semicolon_values <- function(value, numeric = FALSE) {
    value <- .as_ena_scalar(value)
    if (!nzchar(value)) {
        return(if (numeric) numeric() else character())
    }
    values <- trimws(strsplit(value, ";", fixed = TRUE)[[1L]])
    values <- values[nzchar(values)]
    if (!numeric) {
        return(values)
    }
    parsed <- suppressWarnings(as.numeric(values))
    parsed[is.finite(parsed)]
}
.parse_ena_filereport <- function(text, accession) {
    connection <- textConnection(text)
    on.exit(close(connection), add = TRUE)
    raw <- utils::read.delim(connection, check.names = FALSE, stringsAsFactors = FALSE,
        na.strings = c("", "NA"))
    if (!nrow(raw)) {
        stop(sprintf("ENA returned no read-run records for accession %s.", accession),
            call. = FALSE)
    }
    missing <- setdiff(.ena_filereport_fields, colnames(raw))
    if (length(missing)) {
        stop(sprintf("ENA response is missing required fields: %s.", paste(missing,
            collapse = ", ")), call. = FALSE)
    }
    runs <- S4Vectors::DataFrame(run_accession = vapply(raw$run_accession, .as_ena_scalar,
        character(1)), sample_accession = vapply(raw$sample_accession, .as_ena_scalar,
        character(1)), sample_title = vapply(raw$sample_title, .as_ena_scalar,
        character(1)), library_name = vapply(raw$library_name, .as_ena_scalar,
        character(1)), experiment_title = vapply(raw$experiment_title, .as_ena_scalar,
        character(1)), read_count = suppressWarnings(as.numeric(raw$read_count)),
        base_count = suppressWarnings(as.numeric(raw$base_count)), fastq_urls = I(lapply(raw$fastq_ftp,
            .ena_semicolon_values)), fastq_md5 = I(lapply(raw$fastq_md5, .ena_semicolon_values)),
        fastq_bytes = I(lapply(raw$fastq_bytes, .ena_semicolon_values, numeric = TRUE)))
    first <- raw[1L, , drop = FALSE]
    structure(list(accession = accession, bioproject_id = .as_ena_scalar(first$study_accession),
        study_title = .as_ena_scalar(first$study_title), library_strategy = .as_ena_scalar(first$library_strategy),
        library_source = .as_ena_scalar(first$library_source), runs = runs), class = c("selexprep_inspection",
        "list"))
}
.fetch_plan <- function(inspection, overrides = numeric()) {
    if (!inherits(inspection, "selexprep_inspection")) {
        stop("`inspection` must be a selexprep_inspection.", call. = FALSE)
    }
    assignments <- infer_selexprep_rounds(inspection$runs, overrides = overrides)
    rows <- lapply(seq_len(nrow(inspection$runs)), function(index) {
        urls <- inspection$runs$fastq_urls[[index]]
        if (!length(urls)) {
            return(NULL)
        }
        md5 <- inspection$runs$fastq_md5[[index]]
        bytes <- inspection$runs$fastq_bytes[[index]]
        n <- length(urls)
        data.frame(run_accession = rep(inspection$runs$run_accession[[index]],
            n), round_number = rep(assignments$round_number[[index]], n), round_confidence = rep(assignments$confidence[[index]],
            n), round_source = rep(assignments$source_field[[index]], n), round_unambiguous = rep(length(assignments$round_candidates[[index]]) ==
            1L, n), file_name = basename(sub("\\?.*$", "", urls)), url = ifelse(grepl("^https?://",
            urls), urls, paste0("https://", sub("^ftp://", "", urls))), expected_md5 = c(md5,
            rep(NA_character_, n))[seq_len(n)], expected_bytes = c(bytes, rep(NA_real_,
            n))[seq_len(n)], stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) {
        stop("ENA did not provide FASTQ URLs for this accession.", call. = FALSE)
    }
    plan <- S4Vectors::DataFrame(do.call(rbind, rows))
    if (anyDuplicated(plan$file_name)) {
        stop("ENA returned duplicated FASTQ basenames; refusing an ambiguous download plan.",
            call. = FALSE)
    }
    plan
}
#' Inspect public sequencing metadata at ENA
#'
#' Queries ENA's read-run filereport endpoint without downloading sequence data.
#' The returned `selexprep_inspection` retains verbatim library strategy/source
#' metadata and one row per run with its FASTQ URLs, sizes, and MD5 checksums.
#'
#' @param accession Study- or run-level ENA/INSDC accession.
#' @param timeout_seconds Positive network timeout in seconds.
#'
#' @return A `selexprep_inspection` list with a run-level
#'   `S4Vectors::DataFrame` in `$runs`.
#' @export
#' @examples
#' \dontrun{
#' inspection <- inspect_selexprep_accession('PRJDB19098')
#' inspection$runs
#' }
inspect_selexprep_accession <- function(accession, timeout_seconds = 30) {
    if (!is.character(accession) || length(accession) != 1L || is.na(accession) ||
        !nzchar(accession)) {
        stop("`accession` must be one non-empty string.", call. = FALSE)
    }
    if (!is.numeric(timeout_seconds) || length(timeout_seconds) != 1L || is.na(timeout_seconds) ||
        timeout_seconds <= 0) {
        stop("`timeout_seconds` must be one positive number.", call. = FALSE)
    }
    request <- httr2::request(.ena_filereport_url)
    request <- httr2::req_url_query(request, accession = accession, result = "read_run",
        fields = paste(.ena_filereport_fields, collapse = ","), format = "tsv")
    request <- httr2::req_timeout(request, timeout_seconds)
    request <- httr2::req_retry(request, max_tries = 3L)
    request <- httr2::req_error(request, is_error = function(response) FALSE)
    response <- httr2::req_perform(request)
    if (httr2::resp_status(response) >= 400L) {
        stop(sprintf("ENA request failed with HTTP status %d.", httr2::resp_status(response)),
            call. = FALSE)
    }
    .parse_ena_filereport(httr2::resp_body_string(response), accession)
}
#' Fetch FASTQ files from an ENA inspection
#'
#' Creates a deterministic, checksum-aware download plan from an existing
#' inspection or an accession. By default no bytes are transferred; set
#' `dry_run = FALSE` to download to `outdir`. Each completed file is checked
#' against ENA's MD5 value when one is supplied.
#'
#' @param accession_or_inspection A character accession or a
#'   `selexprep_inspection` returned by `inspect_selexprep_accession()`.
#' @param outdir Destination directory for downloaded FASTQ files.
#' @param dry_run Whether to return the plan without downloading.
#' @param overwrite Whether an existing destination file may be replaced.
#' @param timeout_seconds Positive timeout for metadata and FASTQ requests.
#' @param require_assigned_rounds Whether to refuse ambiguous or missing round
#'   assignments before creating a download plan.
#' @param overrides Optional named integer mapping of run accessions to manually
#'   curated round numbers.
#'
#' @return A `selexprep_fetch_result` list containing the run-level `$plan`,
#'   `$downloaded_files`, and `$dry_run`.
#' @export
#' @examples
#' \dontrun{
#' inspection <- inspect_selexprep_accession('PRJDB19098')
#' fetch_selexprep_reads(inspection, 'raw/PRJDB19098')
#' }
fetch_selexprep_reads <- function(accession_or_inspection, outdir, dry_run = TRUE,
    overwrite = FALSE, timeout_seconds = 30, require_assigned_rounds = TRUE, overrides = numeric()) {
    valid_timeout <- is.numeric(timeout_seconds) && length(timeout_seconds) ==
        1L && !is.na(timeout_seconds) && timeout_seconds > 0
    if (!valid_timeout) {
        stop("`timeout_seconds` must be one positive number.", call. = FALSE)
    }
    inspection <- if (is.character(accession_or_inspection)) {
        inspect_selexprep_accession(accession_or_inspection, timeout_seconds = timeout_seconds)
    } else {
        accession_or_inspection
    }
    if (!is.logical(dry_run) || length(dry_run) != 1L || is.na(dry_run) || !is.logical(overwrite) ||
        length(overwrite) != 1L || is.na(overwrite) || !is.logical(require_assigned_rounds) ||
        length(require_assigned_rounds) != 1L || is.na(require_assigned_rounds)) {
        stop("`dry_run`, `overwrite`, and `require_assigned_rounds` must be single non-missing logical values.",
            call. = FALSE)
    }
    if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) || !nzchar(outdir)) {
        stop("`outdir` must be one non-empty path.", call. = FALSE)
    }
    plan <- .fetch_plan(inspection, overrides = overrides)
    if (require_assigned_rounds) {
        unsafe <- is.na(plan$round_number) | !plan$round_unambiguous
        if (any(unsafe)) {
            stop(sprintf("Refusing download: no unambiguous SELEX round assignment for %s. Supply `overrides` after manual review or set `require_assigned_rounds = FALSE` for acquisition only.",
                paste(unique(plan$run_accession[unsafe]), collapse = ", ")), call. = FALSE)
        }
    }
    if (dry_run) {
        return(structure(list(plan = plan, downloaded_files = character(), dry_run = TRUE),
            class = c("selexprep_fetch_result", "list")))
    }
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    downloaded_files <- character(nrow(plan))
    for (index in seq_len(nrow(plan))) {
        destination <- file.path(outdir, plan$file_name[[index]])
        if (file.exists(destination) && !overwrite) {
            stop(sprintf("Refusing to overwrite existing file: %s", destination),
                call. = FALSE)
        }
        temporary <- tempfile(pattern = ".selexprep-download-", tmpdir = outdir)
        request <- httr2::request(plan$url[[index]])
        request <- httr2::req_timeout(request, timeout_seconds)
        request <- httr2::req_retry(request, max_tries = 3L)
        request <- httr2::req_error(request, is_error = function(response) FALSE)
        response <- httr2::req_perform(request, path = temporary)
        if (httr2::resp_status(response) >= 400L) {
            stop(sprintf("Download failed for %s with HTTP status %d.", plan$url[[index]],
                httr2::resp_status(response)), call. = FALSE)
        }
        expected_md5 <- plan$expected_md5[[index]]
        if (!is.na(expected_md5)) {
            observed_md5 <- digest::digest(temporary, algo = "md5", file = TRUE)
            if (!identical(tolower(observed_md5), tolower(expected_md5))) {
                stop(sprintf("MD5 verification failed for %s; temporary file retained at %s.",
                  plan$file_name[[index]], temporary), call. = FALSE)
            }
        }
        if (!file.rename(temporary, destination)) {
            stop(sprintf("Could not move verified download to %s.", destination),
                call. = FALSE)
        }
        downloaded_files[[index]] <- normalizePath(destination, winslash = "/")
    }
    structure(list(plan = plan, downloaded_files = downloaded_files, dry_run = FALSE),
        class = c("selexprep_fetch_result", "list"))
}
