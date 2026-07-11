#' selexprepR result objects
#'
#' The package uses lightweight S3 lists for decision records and reports, while
#' sequence abundance data use the standard Bioconductor `SummarizedExperiment`
#' container. These lists are deliberately explicit and serializable: they can
#' be written to the versioned JSON contracts without loss of meaning.
#'
#' The list elements are stable public contracts:
#'
#' * `selexprep_library_report` stores `primer_5p`, `primer_3p`, their variant
#'   evidence (`variants_5p`, `variants_3p`), `known_adapter_hits`,
#'   `extraction_mode`, `full_insert_recovered`, `read_source`,
#'   `required_action`, `orientation`, random-region length summaries
#'   (`n_length_mode`, `n_length_distribution`, `n_length_confidence`), primer
#'   match and position-consistency rates, the inference fraction and seed,
#'   overall `confidence`, `status`, and an optional `failure_reason`.
#' * `selexprep_extraction` stores `extraction_mode`,
#'   `full_insert_recovered`, primer-stripped `sequences_by_round`, named
#'   `input_reads` and `output_reads` counts, and the pre-trim
#'   forward/reverse/ambiguous `strand_distribution` for each round.
#' * `selexprep_qc` stores a `per_round` metrics `DataFrame`, a pairwise
#'   `round_similarity` `DataFrame`, a `flags` `DataFrame` whose `evidence`
#'   column contains structured lists, and the applied `settings`.
#' * `selexprep_manifest` stores the schema and package/runtime versions,
#'   accession and run identifiers, input/output SHA-256 maps, the embedded
#'   library report and its classification, parameters, stage durations, QC
#'   flags, and sampling seed. Readers accept portable Python v1 and R-native
#'   v2 schemas.
#' * `selexprep_inspection` stores the requested `accession`, BioProject ID,
#'   study title, library strategy/source, and a run-level `runs` `DataFrame`
#'   containing accessions, titles, sizes, FASTQ URLs, and checksums.
#' * `selexprep_fetch_result` stores the resolved run-level download `plan`,
#'   `downloaded_files`, and the `dry_run` status.
#' * `selexprep_read_inputs` stores round-aware R1 pools, optional R2 pools,
#'   the inferred `read_source`, and a file-level provenance `DataFrame`.
#' * `selexprep_demultiplex` extends `selexprep_read_inputs` with unassigned
#'   reads, per-round assignment summaries, and the validated barcode map.
#'
#' S3 lists are used for these small decision and provenance records because
#' their named elements map losslessly to the versioned JSON contracts. The
#' high-dimensional assay remains a standard `SummarizedExperiment` with a
#' sparse `Matrix` assay.
#'
#' @name selexprep-classes
#' @aliases selexprep_library_report
#' @aliases selexprep_extraction
#' @aliases selexprep_manifest
#' @aliases selexprep_inspection
#' @aliases selexprep_fetch_result
#' @aliases selexprep_read_inputs
#' @examples
#' p5 <- 'GGTAATACGACTCACTATAGGG'
#' p3 <- 'CCATGCATGCATGCATGCAT'
#' report <- selexprep_detect(list(round_00 = rep(paste0(p5, 'ACGTACGTACGTACGT', p3), 500)))
#' class(report)
NULL
