.ena_filereport_url <- "https://www.ebi.ac.uk/ena/portal/api/filereport"
.ena_filereport_fields <- c("run_accession", "study_accession", "study_title",
    "library_strategy", "library_source", "library_name", "experiment_title", "sample_title",
    "sample_accession", "read_count", "base_count", "fastq_md5", "fastq_bytes",
    "fastq_ftp")
.ena_sample_xml_url <- "https://www.ebi.ac.uk/ena/browser/api/xml/"
.ena_perform_request <- function(request) httr2::req_perform(request)
.ena_retry_sleep <- function(seconds) Sys.sleep(seconds)
.ena_metadata_request <- function(url, timeout_seconds) {
    request <- httr2::request(url)
    request <- httr2::req_timeout(request, timeout_seconds)
    request <- httr2::req_retry(request, max_tries = 3L)
    httr2::req_error(request, is_error = function(response) FALSE)
}
.ena_filereport_text <- function(accession, timeout_seconds) {
    for (attempt in seq_len(4L)) {
        request <- .ena_metadata_request(.ena_filereport_url, timeout_seconds)
        request <- httr2::req_url_query(request, accession = accession,
            result = "read_run", fields = paste(.ena_filereport_fields, collapse = ","),
            format = "tsv")
        response <- .ena_perform_request(request)
        if (httr2::resp_status(response) >= 400L) {
            stop(sprintf("ENA request failed with HTTP status %d.", httr2::resp_status(response)),
                call. = FALSE)
        }
        text <- httr2::resp_body_string(response)
        if (!grepl("ERROR occurred", text, fixed = TRUE)) return(text)
        if (attempt < 4L) .ena_retry_sleep(2 ^ (attempt - 1L))
    }
    stop(sprintf(
        "ENA filereport did not return records after 4 attempts for %s; %s",
        accession,
        "this is the archive's service, not the data: retry later."
    ), call. = FALSE)
}
.parse_ena_sample_attributes <- function(body, wanted) {
    document <- xml2::read_xml(body)
    samples <- xml2::xml_find_all(document, "//SAMPLE")
    out <- list()
    for (sample in samples) {
        attributes <- character()
        for (attribute in xml2::xml_find_all(sample, ".//SAMPLE_ATTRIBUTE")) {
            tag <- trimws(xml2::xml_text(xml2::xml_find_first(attribute, "TAG")))
            value <- trimws(xml2::xml_text(xml2::xml_find_first(attribute, "VALUE")))
            if (is.na(tag) || !nzchar(tag) || startsWith(tag, "ENA-") || tag %in% names(attributes)) next
            if (is.na(value)) value <- ""
            attributes <- c(attributes, stats::setNames(value, tag))
        }
        identifiers <- c(xml2::xml_attr(sample, "accession"),
            xml2::xml_text(xml2::xml_find_all(sample, ".//PRIMARY_ID | .//EXTERNAL_ID")))
        identifiers <- intersect(trimws(identifiers[!is.na(identifiers)]), wanted)
        for (identifier in identifiers) out[[identifier]] <- attributes
    }
    out
}
.query_ena_sample_attributes <- function(sample_accessions, timeout_seconds) {
    accessions <- sort(unique(trimws(sample_accessions[!is.na(sample_accessions)])))
    accessions <- accessions[nzchar(accessions)]
    out <- list()
    if (!length(accessions)) return(out)
    get <- function(batch) {
        request <- .ena_metadata_request(paste0(.ena_sample_xml_url,
            paste(batch, collapse = ",")), timeout_seconds)
        .ena_perform_request(request)
    }
    collect <- function(response, wanted) {
        status <- httr2::resp_status(response)
        if (status == 404L && length(wanted) == 1L) return(list())
        if (status >= 400L) {
            stop(sprintf("ENA sample request failed with HTTP status %d.", status), call. = FALSE)
        }
        .parse_ena_sample_attributes(httr2::resp_body_raw(response), wanted)
    }
    for (start in seq.int(1L, length(accessions), by = 50L)) {
        batch <- accessions[seq.int(start, min(start + 49L, length(accessions)))]
        response <- get(batch)
        if (httr2::resp_status(response) == 404L && length(batch) > 1L) {
            for (accession in batch) {
                attributes <- collect(get(accession), accession)
                out[names(attributes)] <- attributes
            }
        } else {
            attributes <- collect(response, batch)
            out[names(attributes)] <- attributes
        }
    }
    out
}
.enrich_ena_sample_attributes <- function(inspection, timeout_seconds) {
    attributes <- tryCatch(.query_ena_sample_attributes(inspection$runs$sample_accession,
        timeout_seconds), error = function(error) {
            warning(sprintf(
                "Sample attributes of %s could not be fetched (%s); %s",
                inspection$accession,
                conditionMessage(error),
                "rounds are parsed from titles and library names only."
            ), call. = FALSE)
            list()
        })
    inspection$runs$sample_attributes <- I(lapply(inspection$runs$sample_accession,
        function(sample) {
            if (sample %in% names(attributes)) attributes[[sample]] else character()
        }))
    inspection
}
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
        stop(sprintf("ENA response is missing required fields: %s.", toString(missing)),
            call. = FALSE)
    }
    runs <- S4Vectors::DataFrame(run_accession = vapply(raw$run_accession, .as_ena_scalar,
        character(1)), sample_accession = vapply(raw$sample_accession, .as_ena_scalar,
        character(1)), sample_title = vapply(raw$sample_title, .as_ena_scalar,
        character(1)), library_name = vapply(raw$library_name, .as_ena_scalar,
        character(1)), experiment_title = vapply(raw$experiment_title, .as_ena_scalar,
        character(1)), library_strategy = unname(vapply(raw$library_strategy, .as_ena_scalar,
        character(1))), sample_attributes = I(rep(list(character()), nrow(raw))),
        read_count = suppressWarnings(as.numeric(raw$read_count)),
        base_count = suppressWarnings(as.numeric(raw$base_count)), fastq_urls = I(lapply(raw$fastq_ftp,
            .ena_semicolon_values)), fastq_md5 = I(lapply(raw$fastq_md5, .ena_semicolon_values)),
        fastq_bytes = I(lapply(raw$fastq_bytes, .ena_semicolon_values, numeric = TRUE)))
    first <- raw[1L, , drop = FALSE]
    structure(list(accession = accession, bioproject_id = .as_ena_scalar(first$study_accession),
        study_title = .as_ena_scalar(first$study_title), library_strategy = .as_ena_scalar(first$library_strategy),
        library_source = .as_ena_scalar(first$library_source), runs = runs), class = c("selexprep_inspection",
        "list"))
}
.rounds_with_several_samples <- function(runs, assignments) {
    safe <- .safe_round_assignments(assignments)
    grouped <- split(which(safe), assignments$round_number[safe])
    shared <- lapply(grouped, function(indices) {
        sample <- if ("sample_accession" %in% colnames(runs)) {
            as.character(runs$sample_accession[indices])
        } else {
            rep("", length(indices))
        }
        accession <- as.character(runs$run_accession[indices])
        missing <- is.na(sample) | !nzchar(sample)
        sample[missing] <- accession[missing]
        if (length(unique(sample)) > 1L) accession else NULL
    })
    shared <- Filter(Negate(is.null), shared)
    if (!length(shared)) return(list())
    shared[order(as.numeric(names(shared)))]
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
        row <- data.frame(run_accession = rep(inspection$runs$run_accession[[index]],
            n), round_number = rep(assignments$round_number[[index]], n), round_confidence = rep(assignments$confidence[[index]],
            n), round_source = rep(assignments$source_field[[index]], n), round_unambiguous = rep(length(assignments$round_candidates[[index]]) ==
            1L, n), file_name = basename(sub("\\?.*$", "", urls)), url = ifelse(grepl("^https?://",
            urls), urls, paste0("https://", sub("^ftp://", "", urls))), expected_md5 = c(md5,
            rep(NA_character_, n))[seq_len(n)], expected_bytes = c(bytes, rep(NA_real_,
            n))[seq_len(n)], stringsAsFactors = FALSE)
        row$sample_accession <- rep(if ("sample_accession" %in% colnames(inspection$runs))
            inspection$runs$sample_accession[[index]] else "", n)
        row$library_strategy <- rep(if ("library_strategy" %in% colnames(inspection$runs))
            inspection$runs$library_strategy[[index]] else inspection$library_strategy, n)
        row$sample_attributes <- I(rep(list(if ("sample_attributes" %in% colnames(inspection$runs))
            inspection$runs$sample_attributes[[index]] else character()), n))
        row
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
    S4Vectors::metadata(plan)$rounds_with_several_samples <-
        .rounds_with_several_samples(inspection$runs, assignments)
    plan
}
#' Inspect public sequencing metadata at ENA
#'
#' Queries ENA's read-run filereport endpoint without downloading sequence data.
#' The returned `selexprep_inspection` retains verbatim library strategy/source
#' metadata and one row per run with its FASTQ URLs, sizes, and MD5 checksums.
#' Custom BioSample attributes are fetched from ENA sample XML and retained
#' in a `sample_attributes` list-column; each run keeps its `library_strategy`.
#' If sample XML is unavailable, a warning accompanies fallback to text-field
#' round inference. Filereport service-error bodies are retried before parsing.
#'
#' @param accession Study- or run-level ENA/INSDC accession.
#' @param timeout_seconds Positive network timeout in seconds.
#'
#' @return A `selexprep_inspection` list with a run-level
#'   `S4Vectors::DataFrame` in `$runs`.
#' @export
#' @examples
#' \donttest{
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
    text <- .ena_filereport_text(accession, timeout_seconds)
    inspection <- .parse_ena_filereport(text, accession)
    .enrich_ena_sample_attributes(inspection, timeout_seconds)
}
#' Fetch FASTQ files from an ENA inspection
#'
#' Creates a deterministic, checksum-aware download plan from an existing
#' inspection or an accession. By default no bytes are transferred; set
#' `dry_run = FALSE` to download to `outdir`. Each completed file is checked
#' against ENA's MD5 value when one is supplied.
#' Assigned rounds spanning distinct BioSamples emit a warning: counting
#' merges their runs into one round pool, so review replicates or parallel
#' selections before proceeding. The assignments are retained.
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
#'   `$downloaded_files`, `$dry_run`, and `$rounds_with_several_samples`, a
#'   mapping of assigned rounds to runs from distinct BioSamples. The mapping
#'   is also stored in `S4Vectors::metadata(result$plan)`.
#' @export
#' @examples
#' \donttest{
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
    shared <- S4Vectors::metadata(plan)$rounds_with_several_samples
    if (length(shared)) {
        labels <- vapply(names(shared), function(round) {
            sprintf("round %s: %s", round, paste(shared[[round]], collapse = ", "))
        }, character(1))
        warning(sprintf(paste0("%d round(s) carry runs from more than one BioSample (%s). ",
            "Counting merges every run of a round into one pool: suitable for replicates of one selection, ",
            "but separate selections require separate run selections and round mappings."),
            length(shared), paste(labels, collapse = "; ")), call. = FALSE)
    }
    if (dry_run) {
        return(structure(list(plan = plan, downloaded_files = character(), dry_run = TRUE,
            rounds_with_several_samples = shared),
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
    structure(list(plan = plan, downloaded_files = downloaded_files, dry_run = FALSE,
        rounds_with_several_samples = shared),
        class = c("selexprep_fetch_result", "list"))
}
