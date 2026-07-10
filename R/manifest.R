.manifest_v1_fields <- c(
    "manifest_version", "selexprep_version", "python_version", "cutadapt_version",
    "dnaio_version", "pyarrow_version", "accession", "bioproject_id", "runs",
    "input_sha256", "output_sha256", "library_report", "extraction_mode",
    "read_source", "required_action", "full_insert_recovered", "parameters",
    "runtime_seconds_per_stage", "flags", "sampling_seed"
)

.manifest_v2_fields <- c(
    "manifest_version", "selexprep_version", "r_version", "bioconductor_version",
    "accession", "bioproject_id", "runs", "input_sha256", "output_sha256",
    "library_report", "extraction_mode", "read_source", "required_action",
    "full_insert_recovered", "parameters", "runtime_seconds_per_stage", "flags",
    "sampling_seed"
)

.manifest_fields_for_version <- function(version) {
    switch(version,
        selexprep_manifest_v1 = .manifest_v1_fields,
        selexprep_manifest_v2 = .manifest_v2_fields,
        stop("Unsupported selexprep manifest version.", call. = FALSE)
    )
}

.as_optional_string <- function(value, name) {
    if (is.null(value)) {
        return(NULL)
    }
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
        stop(sprintf("`%s` must be NULL or one non-missing string.", name), call. = FALSE)
    }
    unname(value)
}

.as_string <- function(value, name) {
    value <- .as_optional_string(value, name)
    if (is.null(value)) {
        stop(sprintf("`%s` must be one non-missing string.", name), call. = FALSE)
    }
    value
}

.as_named_map <- function(value, name, mode = c("character", "numeric")) {
    mode <- match.arg(mode)
    if (is.null(value)) {
        stop(sprintf("`%s` must be a named map, not NULL.", name), call. = FALSE)
    }
    if (is.list(value)) {
        value <- unlist(value, use.names = TRUE)
    }
    if (!length(value)) {
        return(if (identical(mode, "character")) character() else numeric())
    }
    if (is.null(names(value)) || any(!nzchar(names(value))) || anyDuplicated(names(value))) {
        stop(sprintf("`%s` must have unique, non-empty names.", name), call. = FALSE)
    }
    if (identical(mode, "character")) {
        if (!is.character(value) || anyNA(value)) {
            stop(sprintf("`%s` must contain non-missing strings.", name), call. = FALSE)
        }
        return(stats::setNames(unname(value), names(value)))
    }
    if (!is.numeric(value) || anyNA(value) || any(!is.finite(value))) {
        stop(sprintf("`%s` must contain finite numeric values.", name), call. = FALSE)
    }
    stats::setNames(as.numeric(value), names(value))
}

.as_character_vector <- function(value, name) {
    if (is.list(value)) {
        value <- unlist(value, use.names = FALSE)
    }
    if (is.null(value) || !length(value)) {
        return(character())
    }
    if (!is.character(value) || anyNA(value)) {
        stop(sprintf("`%s` must be a character vector without missing values.", name), call. = FALSE)
    }
    unname(value)
}

.as_sorted_json_map <- function(value) {
    if (!length(value)) {
        return(list())
    }
    as.list(value[order(names(value), method = "radix")])
}

.library_report_from_payload <- function(payload) {
    if (!is.list(payload) || !identical(names(payload), .library_report_fields)) {
        stop("`library_report` does not match the stable LibraryReport schema.", call. = FALSE)
    }
    payload$known_adapter_hits <- .as_named_map(
        payload$known_adapter_hits,
        "library_report$known_adapter_hits",
        mode = "numeric"
    )
    payload$n_length_distribution <- .as_named_map(
        payload$n_length_distribution,
        "library_report$n_length_distribution",
        mode = "numeric"
    )
    .new_library_report(payload)
}

.selexprep_version <- function() {
    version <- tryCatch(
        utils::packageDescription("selexprep", fields = "Version"),
        error = function(...) NULL
    )
    if (is.null(version) || !length(version) || is.na(version)) {
        return("0.99.0")
    }
    as.character(version)
}

.bioconductor_version <- function() {
    version <- tryCatch(utils::packageVersion("BiocVersion"), error = function(...) NULL)
    if (is.null(version)) NULL else as.character(version)
}

.new_selexprep_manifest <- function(fields) {
    if (!is.list(fields) || is.null(fields$manifest_version)) {
        stop("A manifest must be a named list with `manifest_version`.", call. = FALSE)
    }
    version <- .as_string(fields$manifest_version, "manifest_version")
    expected_fields <- .manifest_fields_for_version(version)
    if (!identical(names(fields), expected_fields)) {
        stop("The manifest does not contain the stable field set for its version.", call. = FALSE)
    }
    fields$manifest_version <- version
    fields$selexprep_version <- .as_string(fields$selexprep_version, "selexprep_version")
    if (identical(version, "selexprep_manifest_v1")) {
        for (name in c("python_version", "cutadapt_version", "dnaio_version", "pyarrow_version")) {
            fields[[name]] <- .as_string(fields[[name]], name)
        }
    } else {
        fields$r_version <- .as_string(fields$r_version, "r_version")
        fields["bioconductor_version"] <- list(.as_optional_string(
            fields$bioconductor_version,
            "bioconductor_version"
        ))
    }
    fields["accession"] <- list(.as_optional_string(fields$accession, "accession"))
    fields["bioproject_id"] <- list(.as_optional_string(fields$bioproject_id, "bioproject_id"))
    fields$runs <- .as_character_vector(fields$runs, "runs")
    fields$input_sha256 <- .as_named_map(fields$input_sha256, "input_sha256")
    fields$output_sha256 <- .as_named_map(fields$output_sha256, "output_sha256")
    fields$library_report <- .library_report_from_payload(fields$library_report)
    for (name in c("extraction_mode", "read_source", "required_action")) {
        fields[[name]] <- .as_string(fields[[name]], name)
    }
    if (!identical(fields$extraction_mode, fields$library_report$extraction_mode) ||
        !identical(fields$read_source, fields$library_report$read_source) ||
        !identical(fields$required_action, fields$library_report$required_action)) {
        stop("Manifest classification fields must agree with `library_report`.", call. = FALSE)
    }
    if (!is.logical(fields$full_insert_recovered) || length(fields$full_insert_recovered) != 1L ||
        is.na(fields$full_insert_recovered) ||
        !identical(fields$full_insert_recovered, fields$library_report$full_insert_recovered)) {
        stop("`full_insert_recovered` must agree with `library_report`.", call. = FALSE)
    }
    fields$parameters <- .as_named_map(fields$parameters, "parameters")
    fields$runtime_seconds_per_stage <- .as_named_map(
        fields$runtime_seconds_per_stage,
        "runtime_seconds_per_stage",
        mode = "numeric"
    )
    fields$flags <- .as_character_vector(fields$flags, "flags")
    if (!is.numeric(fields$sampling_seed) || length(fields$sampling_seed) != 1L ||
        is.na(fields$sampling_seed) || !identical(as.numeric(fields$sampling_seed),
            as.numeric(fields$library_report$sampling_seed))) {
        stop("`sampling_seed` must agree with `library_report`.", call. = FALSE)
    }
    fields$sampling_seed <- as.numeric(fields$sampling_seed)

    structure(fields, class = c("selexprep_manifest", "list"))
}

.manifest_payload <- function(manifest) {
    if (!inherits(manifest, "selexprep_manifest")) {
        stop("`manifest` must be a selexprep_manifest.", call. = FALSE)
    }
    payload <- unclass(manifest)
    payload$input_sha256 <- .as_sorted_json_map(payload$input_sha256)
    payload$output_sha256 <- .as_sorted_json_map(payload$output_sha256)
    payload$parameters <- .as_sorted_json_map(payload$parameters)
    payload$runtime_seconds_per_stage <- .as_sorted_json_map(payload$runtime_seconds_per_stage)
    payload$library_report <- .library_report_payload(payload$library_report)
    payload
}

.hash_paths <- function(paths, root = NULL) {
    paths <- as.character(paths)
    if (!length(paths)) {
        return(character())
    }
    expanded <- normalizePath(paths, winslash = "/", mustWork = TRUE)
    if (!is.null(root)) {
        root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    }
    keys <- vapply(expanded, function(path) {
        if (!is.null(root) && startsWith(path, paste0(root, "/"))) {
            substr(path, nchar(root) + 2L, nchar(path))
        } else {
            basename(path)
        }
    }, character(1))
    if (anyDuplicated(keys)) {
        stop("Paths produce duplicated manifest hash keys.", call. = FALSE)
    }
    hashes <- vapply(expanded, digest::digest, character(1), algo = "sha256", file = TRUE)
    stats::setNames(unname(hashes), keys)
}

#' Build an R-native reproducibility manifest
#'
#' Builds a versioned `selexprep_manifest_v2` object. The reader also accepts
#' the historical Python `selexprep_manifest_v1` schema, so existing analysis
#' records remain inspectable. Hashes are computed for every supplied existing
#' path and are keyed relative to `input_root` or `output_root` when supplied.
#'
#' @param library_report A `selexprep_library_report`.
#' @param input_paths Existing input file paths to hash.
#' @param output_paths Existing output file paths to hash.
#' @param accession Optional public accession.
#' @param bioproject_id Optional BioProject accession.
#' @param runs Run accessions.
#' @param parameters Named character parameters used for the run.
#' @param runtime_seconds_per_stage Named numeric stage durations.
#' @param flags Character QC flag names.
#' @param input_root Optional root used to make input hash keys relative.
#' @param output_root Optional root used to make output hash keys relative.
#'
#' @return A validated `selexprep_manifest`.
#' @export
#' @examples
#' p5 <- "GGTAATACGACTCACTATAGGG"
#' p3 <- "CCATGCATGCATGCATGCAT"
#' report <- selexprep_detect(list(round_00 = rep(paste0(p5, "ACGTACGTACGTACGT", p3), 500)))
#' build_selexprep_manifest(report)
build_selexprep_manifest <- function(library_report, input_paths = character(),
    output_paths = character(), accession = NULL, bioproject_id = NULL,
    runs = character(), parameters = character(), runtime_seconds_per_stage = numeric(),
    flags = character(), input_root = NULL, output_root = NULL) {
    if (!inherits(library_report, "selexprep_library_report")) {
        stop("`library_report` must be a selexprep_library_report.", call. = FALSE)
    }
    .new_selexprep_manifest(list(
        manifest_version = "selexprep_manifest_v2",
        selexprep_version = .selexprep_version(),
        r_version = as.character(getRversion()),
        bioconductor_version = .bioconductor_version(),
        accession = accession,
        bioproject_id = bioproject_id,
        runs = runs,
        input_sha256 = .hash_paths(input_paths, root = input_root),
        output_sha256 = .hash_paths(output_paths, root = output_root),
        library_report = .library_report_payload(library_report),
        extraction_mode = library_report$extraction_mode,
        read_source = library_report$read_source,
        required_action = library_report$required_action,
        full_insert_recovered = library_report$full_insert_recovered,
        parameters = parameters,
        runtime_seconds_per_stage = runtime_seconds_per_stage,
        flags = flags,
        sampling_seed = library_report$sampling_seed
    ))
}

#' Write a selexprep manifest as deterministic JSON
#'
#' @param manifest A `selexprep_manifest`.
#' @param path Output JSON path.
#'
#' @return The SHA-256 digest of the emitted UTF-8 JSON, invisibly.
#' @export
#' @examples
#' p5 <- "GGTAATACGACTCACTATAGGG"
#' p3 <- "CCATGCATGCATGCATGCAT"
#' report <- selexprep_detect(list(round_00 = rep(paste0(p5, "ACGTACGTACGTACGT", p3), 500)))
#' path <- tempfile(fileext = ".json")
#' write_selexprep_manifest(build_selexprep_manifest(report), path)
write_selexprep_manifest <- function(manifest, path) {
    payload <- .manifest_payload(manifest)
    json <- jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA)
    for (name in c("input_sha256", "output_sha256", "parameters", "runtime_seconds_per_stage")) {
        json <- sub(sprintf('"%s": []', name), sprintf('"%s": {}', name), json, fixed = TRUE)
    }
    text <- paste0(json, "\n")
    writeLines(json, con = path, useBytes = TRUE)
    invisible(digest::digest(text, algo = "sha256", serialize = FALSE))
}

#' Read a selexprep manifest JSON file
#'
#' Reads and validates both the portable Python `selexprep_manifest_v1` format
#' and the R-native `selexprep_manifest_v2` format.
#'
#' @param path Path to a selexprep manifest JSON file.
#'
#' @return A validated `selexprep_manifest`.
#' @export
#' @examples
#' p5 <- "GGTAATACGACTCACTATAGGG"
#' p3 <- "CCATGCATGCATGCATGCAT"
#' report <- selexprep_detect(list(round_00 = rep(paste0(p5, "ACGTACGTACGTACGT", p3), 500)))
#' path <- tempfile(fileext = ".json")
#' write_selexprep_manifest(build_selexprep_manifest(report), path)
#' read_selexprep_manifest(path)
read_selexprep_manifest <- function(path) {
    payload <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    .new_selexprep_manifest(payload)
}
