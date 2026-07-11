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
#' selexprep_catalog('FGF-9')
selexprep_catalog <- function(query = NULL) {
    valid_query <- is.character(query) && length(query) == 1L && !is.na(query) &&
        nzchar(query)
    if (!is.null(query) && !valid_query) {
        stop("`query` must be NULL or one non-empty string.", call. = FALSE)
    }
    catalog <- S4Vectors::DataFrame(.load_selexprep_public_catalog())
    if (!is.null(query)) {
        searchable <- lapply(catalog[, c("bioproject_id", "study_title", "target")],
            function(value) {
                value <- as.character(value)
                value[is.na(value)] <- ""
                value
            })
        text <- do.call(paste, searchable)
        keep <- grepl(tolower(query), tolower(text), fixed = TRUE)
        catalog <- catalog[keep, , drop = FALSE]
    }
    catalog
}
.catalog_cache <- new.env(parent = emptyenv())
.load_selexprep_public_catalog <- function() {
    object_name <- "selexprep_public_catalog"
    if (!exists(object_name, envir = .catalog_cache, inherits = FALSE)) {
        package <- utils::packageName(environment())
        utils::data(list = object_name, package = package, envir = .catalog_cache)
    }
    if (!exists(object_name, envir = .catalog_cache, inherits = FALSE)) {
        stop("The bundled selexprepR catalog is unavailable.", call. = FALSE)
    }
    get(object_name, envir = .catalog_cache, inherits = FALSE)
}
