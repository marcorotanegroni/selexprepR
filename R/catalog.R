#' Query the bundled public HT-SELEX study catalog
#'
#' Returns the curated catalog snapshot distributed with the package. The
#' optional text query is matched case-insensitively against the BioProject ID,
#' study title, and target, making it useful for discovery while leaving raw
#' metadata values intact.
#'
#' @param query Optional non-empty search string.
#'
#' @return An `S4Vectors::DataFrame` of curated study metadata.
#' @export
#' @examples
#' selexprep_catalog("FGF-9")
selexprep_catalog <- function(query = NULL) {
    if (!is.null(query) &&
        (!is.character(query) || length(query) != 1L || is.na(query) || !nzchar(query))) {
        stop("`query` must be NULL or one non-empty string.", call. = FALSE)
    }
    path <- system.file("extdata", "selexprep_catalog.csv", package = "selexprep")
    if (!nzchar(path)) {
        stop("The bundled selexprep catalog is unavailable.", call. = FALSE)
    }
    catalog <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
        na.strings = ""
    )
    if (!is.null(query)) {
        searchable <- catalog[, c("bioproject_id", "study_title", "target"), drop = FALSE]
        text <- apply(searchable, 1L, paste, collapse = " ")
        catalog <- catalog[grepl(tolower(query), tolower(text), fixed = TRUE), , drop = FALSE]
    }
    S4Vectors::DataFrame(catalog, check.names = FALSE)
}
