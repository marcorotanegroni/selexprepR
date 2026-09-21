test_that("bundled catalog is available offline and can be queried", {
    catalog <- selexprep_catalog()
    fgf <- selexprep_catalog("FGF-9")

    expect_s4_class(catalog, "DFrame")
    expect_identical(dim(catalog), c(240L, 19L))
    expect_false(anyDuplicated(catalog$bioproject_id) > 0L)
    expect_identical(as.character(fgf$bioproject_id[[1L]]), "PRJDB19098")
    statuses <- as.matrix(as.data.frame(catalog[grep("_curation$", names(catalog))]))
    expect_identical(sum(statuses == "adjudicated"), 47L)
    expect_false(any(statuses == "discordant"))
})

test_that("package dataset and accessor expose the reviewed snapshot", {
    dataset_environment <- new.env(parent = emptyenv())
    utils::data(
        list = "selexprep_public_catalog",
        package = "selexprepR",
        envir = dataset_environment
    )
    dataset <- get(
        "selexprep_public_catalog",
        envir = dataset_environment,
        inherits = FALSE
    )
    provenance <- S4Vectors::metadata(dataset)

    expect_s4_class(dataset, "DFrame")
    expect_identical(dataset, selexprep_catalog())
    expect_identical(
        provenance$snapshot_version,
        "v0.3.2-dual-extraction-adjudicated-en-2026-09-04"
    )
    expect_identical(
        provenance$source_sha256,
        paste0(
            "e67b7b77f9d60b6ae6a686099c2a524ae59f50fbaa446aece7ea607f",
            "dc007fef"
        )
    )
})
