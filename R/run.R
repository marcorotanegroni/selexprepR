#' Run the in-memory SELEX preprocessing workflow
#'
#' Composes primer detection, extraction, sequence counting and quality control
#' into a single Bioconductor result. The primary output is a
#' `SummarizedExperiment` with one column per round and a sparse `counts`
#' assay. The LibraryReport, extraction result and QC report are retained in
#' `metadata()` for reproducibility and inspection.
#'
#' Split-primer paired-end inputs return a valid zero-row experiment because
#' full inserts cannot be counted without read merging. Their separate trimmed
#' sides remain available in `metadata(result)$extraction`.
#'
#' @param sequences_by_round A named list of R1 pools, as accepted by
#'   `selexprep_detect()`, or a `selexprep_read_inputs` object.
#' @param read_source Input read layout.
#' @param paired_mate_streams Optional named R2 pools.
#' @param accession Optional public-study accession stored in metadata.
#' @param sampling_seed Seed forwarded to `selexprep_detect()`.
#' @param max_reads_per_round Optional detection subsampling cap.
#' @param low_total_reads QC threshold for per-round read depth.
#' @param primer_5p Optional manually reviewed 5-prime primer override.
#' @param primer_3p Optional manually reviewed 3-prime primer override.
#' @param extraction_mode Optional manually reviewed extraction-mode override.
#'
#' @return A `SummarizedExperiment` with the LibraryReport, extraction and QC
#'   result in its metadata.
#' @export
#' @examples
#' p5 <- 'GGTAATACGACTCACTATAGGG'
#' p3 <- 'CCATGCATGCATGCATGCAT'
#' pools <- list(round_00 = rep(paste0(p5, 'ACGTACGTACGTACGT', p3), 500))
#' run_selexprep(pools, low_total_reads = 0)
run_selexprep <- function(sequences_by_round, read_source = c("R1", "R2", "R1_AND_R2",
    "INTERLEAVED", "UNKNOWN"), paired_mate_streams = NULL, accession = NULL, sampling_seed = 42L,
    max_reads_per_round = NULL, low_total_reads = 10000L,
    primer_5p = NULL, primer_3p = NULL, extraction_mode = NULL) {
    input_bundle <- NULL
    if (inherits(sequences_by_round, "selexprep_read_inputs")) {
        if (!missing(read_source) || !is.null(paired_mate_streams)) {
            stop(
                "Sequence layout is already defined by the input bundle.",
                call. = FALSE
            )
        }
        input_bundle <- sequences_by_round
        paired_mate_streams <- input_bundle$paired_mate_streams
        read_source <- input_bundle$read_source
        sequences_by_round <- input_bundle$sequences_by_round
    } else {
        read_source <- match.arg(read_source)
    }
    detected <- system.time(report <- selexprep_detect(sequences_by_round = sequences_by_round,
        read_source = read_source, paired_mate_streams = paired_mate_streams, sampling_seed = sampling_seed,
        max_reads_per_round = max_reads_per_round))
    manual_override <- !is.null(primer_5p) || !is.null(primer_3p) ||
        !is.null(extraction_mode)
    report <- .apply_primer_overrides(
        report,
        primer_5p = primer_5p,
        primer_3p = primer_3p,
        extraction_mode = extraction_mode
    )
    extracted <- system.time(extraction <- selexprep_extract(sequences_by_round = sequences_by_round,
        library_report = report, paired_mate_streams = paired_mate_streams))
    counted <- system.time(count_tables <- if (identical(report$extraction_mode,
        "PAIRED_END_SPLIT_PRIMERS")) {
        round_names <- names(extraction$sequences_by_round)
        count_tables <- vector("list", length(round_names))
        for (index in seq_along(count_tables)) {
            count_tables[[index]] <- selexprep_count(character())
        }
        stats::setNames(count_tables, round_names)
    } else {
        lapply(extraction$sequences_by_round, selexprep_count)
    })
    experiment <- .as_selexprep_experiment(count_tables, accession = accession)
    metadata <- S4Vectors::metadata(experiment)
    metadata$library_report <- report
    metadata$extraction <- extraction
    S4Vectors::metadata(experiment) <- metadata
    quality_controlled <- system.time(qc <- if (identical(report$extraction_mode,
        "PAIRED_END_SPLIT_PRIMERS")) {
        NULL
    } else {
        selexprep_qc(experiment, library_report = report, low_total_reads = low_total_reads)
    })
    qc_flags <- if (is.null(qc)) {
        if (identical(report$required_action, "READ_MERGING_RECOMMENDED")) {
            "requires_read_merging_for_full_insert"
        } else {
            character()
        }
    } else {
        as.character(qc$flags$name)
    }
    manifest <- build_selexprep_manifest(
        report,
        accession = accession,
        parameters = c(
            read_source = read_source,
            sampling_seed = as.character(sampling_seed),
            max_reads_per_round = if (is.null(max_reads_per_round)) {
                "all"
            } else {
                as.character(max_reads_per_round)
            },
            low_total_reads = as.character(low_total_reads),
            manual_primer_override = if (manual_override) "true" else "false"
        ),
        runtime_seconds_per_stage = c(
            detect = unname(detected[["elapsed"]]),
            extract = unname(extracted[["elapsed"]]),
            count = unname(counted[["elapsed"]]),
            qc = unname(quality_controlled[["elapsed"]])
        ),
        flags = qc_flags
    )
    metadata <- S4Vectors::metadata(experiment)
    metadata$qc <- qc
    metadata$manifest <- manifest
    metadata$counting_skipped <- identical(report$extraction_mode, "PAIRED_END_SPLIT_PRIMERS")
    if (!is.null(input_bundle) &&
        inherits(input_bundle, "selexprep_demultiplex")) {
        metadata$demultiplex <- list(
            summary = input_bundle$summary,
            barcodes = input_bundle$barcodes,
            max_mismatches = input_bundle$max_mismatches,
            trim_barcode = input_bundle$trim_barcode,
            unassigned_reads = length(input_bundle$unassigned$r1)
        )
    }
    S4Vectors::metadata(experiment) <- metadata
    experiment
}
