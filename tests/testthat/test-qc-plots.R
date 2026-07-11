.make_qc_plot_fixture <- function() {
    as_experiment <- get(
        ".as_selexprep_experiment",
        envir = environment(selexprep_qc)
    )
    experiment <- as_experiment(list(
        round_00 = selexprep_count(c(
            rep("AAAA", 5),
            rep("CCCCCC", 2)
        )),
        round_01 = selexprep_count(c(
            rep("AAAA", 4),
            rep("GGGGG", 3),
            "TTTTTT"
        ))
    ))
    report <- read_library_report(test_path(
        "fixtures",
        "library-report-both-primers.json"
    ))
    extraction <- structure(
        list(
            extraction_mode = "BOTH_PRIMERS_SINGLE_READ",
            full_insert_recovered = TRUE,
            sequences_by_round = list(),
            input_reads = c(round_00 = 10, round_01 = 12),
            output_reads = c(round_00 = 7, round_01 = 8),
            strand_distribution = list()
        ),
        class = c("selexprep_extraction", "list")
    )
    metadata <- S4Vectors::metadata(experiment)
    metadata$extraction <- extraction
    metadata$library_report <- report
    S4Vectors::metadata(experiment) <- metadata

    list(
        experiment = experiment,
        qc = selexprep_qc(
            experiment,
            library_report = report,
            low_total_reads = 0
        )
    )
}

.render_qc_pdf <- function(path, qc, experiment = NULL, ...) {
    plot_method <- get(
        "plot.selexprep_qc",
        envir = environment(selexprep_qc)
    )
    grDevices::pdf(path, width = 11, height = 7)
    on.exit(grDevices::dev.off(), add = TRUE)
    before <- graphics::par("mfrow")
    result <- withVisible(plot_method(qc, y = experiment, ...))
    after <- graphics::par("mfrow")
    list(result = result, before = before, after = after)
}

.pdf_header <- function(path) {
    connection <- file(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    rawToChar(readBin(connection, what = "raw", n = 4L))
}

test_that("QC plot renders the complete diagnostic surface", {
    fixture <- .make_qc_plot_fixture()
    path <- tempfile(fileext = ".pdf")

    rendered <- .render_qc_pdf(
        path,
        fixture$qc,
        fixture$experiment
    )

    expect_false(rendered$result$visible)
    expect_identical(rendered$result$value, fixture$qc)
    expect_identical(rendered$after, rendered$before)
    expect_identical(.pdf_header(path), "%PDF")
    expect_gt(unname(file.info(path)$size), 1000)
})

test_that("QC plot degrades explicitly without an experiment", {
    fixture <- .make_qc_plot_fixture()
    path <- tempfile(fileext = ".pdf")

    expect_silent(.render_qc_pdf(
        path,
        fixture$qc,
        which = c("retention", "primer_match", "length")
    ))

    expect_identical(.pdf_header(path), "%PDF")
    expect_gt(unname(file.info(path)$size), 500)
})

test_that("QC plot validates its inputs and diagnostic names", {
    fixture <- .make_qc_plot_fixture()
    plot_method <- get(
        "plot.selexprep_qc",
        envir = environment(selexprep_qc)
    )

    expect_error(
        plot_method(list()),
        "valid selexprep_qc",
        fixed = TRUE
    )
    expect_error(
        plot_method(fixture$qc, y = data.frame()),
        "SummarizedExperiment",
        fixed = TRUE
    )
    expect_error(
        plot_method(fixture$qc, which = "unknown"),
        "should be one of"
    )
})
