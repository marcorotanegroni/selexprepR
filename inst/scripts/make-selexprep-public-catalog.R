catalog_path <- file.path("inst", "extdata", "selexprep_catalog.csv")
description_path <- "DESCRIPTION"

if (!file.exists(description_path) || !file.exists(catalog_path)) {
    stop("Run this script from the selexprepR package source root.")
}

expected_sha256 <- paste0(
    "9e11bc3868816e769da3ae50c159bb14e679e790be6de106b0ede3756",
    "caa1cb9"
)
expected_columns <- c(
    "bioproject_id", "source", "study_title",
    "study_type", "study_type_curation",
    "target", "target_curation",
    "target_class", "target_class_curation",
    "chemistry", "chemistry_curation",
    "n_random", "n_random_curation",
    "n_rounds", "n_rounds_curation",
    "selection_format", "selection_format_curation",
    "counter_selection", "counter_selection_curation"
)
expected_statuses <- c(
    "concordant", "discordant", "not_stated", "verified",
    "single_source:claude", "single_source:codex"
)

source_sha256 <- digest::digest(catalog_path, algo = "sha256", file = TRUE)
if (!identical(source_sha256, expected_sha256)) {
    stop("Catalog SHA-256 does not match the reviewed source snapshot.")
}

catalog <- utils::read.csv(
    catalog_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "",
    fileEncoding = "UTF-8"
)
catalog[] <- lapply(catalog, function(value) {
    if (is.character(value)) {
        Encoding(value[!is.na(value)]) <- "UTF-8"
    }
    value
})
if (!identical(names(catalog), expected_columns)) {
    stop("Catalog columns differ from the reviewed schema.")
}
if (nrow(catalog) != 240L || anyDuplicated(catalog$bioproject_id)) {
    stop("Catalog rows are incomplete or contain duplicate identifiers.")
}

status_columns <- grep("_curation$", names(catalog), value = TRUE)
observed_statuses <- sort(unique(unlist(catalog[status_columns])))
if (!identical(observed_statuses, sort(expected_statuses))) {
    stop("Catalog contains an unexpected curation status.")
}

selexprep_public_catalog <- S4Vectors::DataFrame(
    catalog,
    check.names = FALSE
)
S4Vectors::metadata(selexprep_public_catalog) <- list(
    snapshot_version = "v0.2.1-dual-extraction-2026-07-03",
    snapshot_date = "2026-07-03",
    source_repository = "https://github.com/marcorotanegroni/selexprep",
    source_commit = "b6792637ec9bed78f131e8f63e12e155f733e9ae",
    source_file = paste0(
        "src/selexprep/catalog/data/curated_metadata.csv"
    ),
    source_sha256 = source_sha256,
    source_url = paste0(
        "https://raw.githubusercontent.com/marcorotanegroni/selexprep/",
        "b6792637ec9bed78f131e8f63e12e155f733e9ae/",
        "src/selexprep/catalog/data/curated_metadata.csv"
    ),
    provenance_url = paste0(
        "https://raw.githubusercontent.com/marcorotanegroni/selexprep/",
        "b6792637ec9bed78f131e8f63e12e155f733e9ae/",
        "src/selexprep/catalog/data/curated_metadata.json"
    ),
    curation_method = paste(
        "Independent Claude and Codex/GPT extraction, followed by",
        "status-preserving reconciliation."
    ),
    terms_notice = paste(
        "The package code is MIT-licensed. Public record metadata remains",
        "subject to the terms of its originating repository or publication."
    )
)

dir.create("data", showWarnings = FALSE)
save(
    selexprep_public_catalog,
    file = file.path("data", "selexprep_public_catalog.rda"),
    compress = "xz",
    version = 3L
)
