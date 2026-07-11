#' Curated public HT-SELEX catalog
#'
#' A frozen, offline snapshot of 240 public HT-SELEX deposits. Eight
#' experimental fields were extracted independently by two language models,
#' reconciled, and retained with a per-field curation status. The flat dataset
#' intentionally omits evidence quotations; immutable links to the canonical
#' provenance-rich JSON and its generation materials are stored in
#' `S4Vectors::metadata(selexprep_public_catalog)`.
#'
#' Missing values mean that a field was not stated in the curated sources. A
#' `discordant` status means that both extracted values are retained in the
#' corresponding value column, separated by `' || '`.
#'
#' @format An `S4Vectors::DataFrame` with 240 rows and 19 columns:
#' \describe{
#'   \item{`bioproject_id`}{ENA BioProject accession or a platform-prefixed
#'     Figshare or Zenodo deposit identifier.}
#'   \item{`source`}{`ena`, or the platform-prefixed Figshare or Zenodo
#'     deposit identifier used as the source key.}
#'   \item{`study_title`}{Title recorded for the public study or deposit.}
#'   \item{`study_type`, `study_type_curation`}{Curated experiment type and
#'     its reconciliation status.}
#'   \item{`target`, `target_curation`}{Selection target and its
#'     reconciliation status.}
#'   \item{`target_class`, `target_class_curation`}{Target class and its
#'     reconciliation status.}
#'   \item{`chemistry`, `chemistry_curation`}{DNA or RNA library chemistry
#'     and its reconciliation status.}
#'   \item{`n_random`, `n_random_curation`}{Reported random-region length and
#'     its reconciliation status.}
#'   \item{`n_rounds`, `n_rounds_curation`}{Reported selection-round count and
#'     its reconciliation status.}
#'   \item{`selection_format`, `selection_format_curation`}{Reported selection
#'     format and its reconciliation status.}
#'   \item{`counter_selection`, `counter_selection_curation`}{Reported
#'     counter-selection conditions and their reconciliation status.}
#' }
#'
#' Curation statuses are `concordant`, `discordant`, `not_stated`, `verified`,
#' `single_source:claude`, or `single_source:codex`.
#'
#' @source The snapshot was migrated without modification from the original
#'   selexprep catalog at commit
#'   `b6792637ec9bed78f131e8f63e12e155f733e9ae`. See the object metadata and
#'   `inst/CATALOG_PROVENANCE.md` for immutable source links and terms.
#' @usage data(selexprep_public_catalog)
#' @keywords datasets
"selexprep_public_catalog"
