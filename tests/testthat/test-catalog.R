test_that("bundled catalog is available offline and can be queried", {
    catalog <- selexprep_catalog()
    fgf <- selexprep_catalog("FGF-9")

    expect_s4_class(catalog, "DFrame")
    expect_gt(nrow(catalog), 200L)
    expect_true("bioproject_id" %in% colnames(catalog))
    expect_identical(as.character(fgf$bioproject_id[[1L]]), "PRJDB19098")
})
