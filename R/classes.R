#' selexprep result objects
#'
#' The package uses lightweight S3 lists for decision records and reports, while
#' sequence abundance data use the standard Bioconductor `SummarizedExperiment`
#' container. These lists are deliberately explicit and serializable: they can
#' be written to the versioned JSON contracts without loss of meaning.
#'
#' `selexprep_library_report` contains the inferred primers, evidence metrics,
#' extraction mode, orientation, and safe-failure reason when inference cannot
#' support extraction. `selexprep_extraction` contains the primer-stripped
#' sequences by round together with input/output read counts and the pre-trim
#' forward/reverse/ambiguous strand distribution for every round.
#'
#' `selexprep_qc` has `per_round` and `flags` `S4Vectors::DataFrame` elements.
#' `selexprep_manifest` stores package/runtime provenance, hashes, the embedded
#' library report, parameters, durations, and QC flags. The reader accepts the
#' portable Python v1 and R-native v2 manifest schemas.
#'
#' `selexprep_inspection` represents ENA metadata and has a run-level `runs`
#' `DataFrame`; `selexprep_fetch_result` contains the corresponding download
#' plan, downloaded paths, and dry-run status.
#'
#' @name selexprep-classes
#' @aliases selexprep_library_report
#' @aliases selexprep_extraction
#' @aliases selexprep_qc
#' @aliases selexprep_manifest
#' @aliases selexprep_inspection
#' @aliases selexprep_fetch_result
#' @examples
#' p5 <- "GGTAATACGACTCACTATAGGG"
#' p3 <- "CCATGCATGCATGCATGCAT"
#' report <- selexprep_detect(list(round_00 = rep(paste0(p5, "ACGTACGTACGTACGT", p3), 500)))
#' class(report)
NULL
