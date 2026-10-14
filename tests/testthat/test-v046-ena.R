.v046_ena_response <- function(text = "", status = 200L) {
    httr2::response(status_code = status, body = charToRaw(text),
        headers = list(`content-type` = "text/plain; charset=UTF-8"))
}
.v046_sample_xml <- paste0(
    '<SAMPLE_SET><SAMPLE accession="SAMEA13183932">',
    '<IDENTIFIERS><PRIMARY_ID>SAMEA13183932</PRIMARY_ID>',
    '<EXTERNAL_ID>SRS000001</EXTERNAL_ID></IDENTIFIERS><SAMPLE_ATTRIBUTES>',
    '<SAMPLE_ATTRIBUTE><TAG>ENA-FIRST-PUBLIC</TAG><VALUE>2022-03-01</VALUE></SAMPLE_ATTRIBUTE>',
    '<SAMPLE_ATTRIBUTE><TAG>construct</TAG><VALUE>ZFC1 &amp; Sall4</VALUE></SAMPLE_ATTRIBUTE>',
    '<SAMPLE_ATTRIBUTE><TAG>selex cycle</TAG><VALUE>1</VALUE></SAMPLE_ATTRIBUTE>',
    '<SAMPLE_ATTRIBUTE><TAG>selex cycle</TAG><VALUE>9</VALUE></SAMPLE_ATTRIBUTE>',
    '</SAMPLE_ATTRIBUTES></SAMPLE></SAMPLE_SET>'
)
.v046_filereport <- function() {
    values <- c("ERR8972349", "PRJEB51212", "Example HT-SELEX", "SELEX", "OTHER",
        "RV04_s", "experiment", "RV04", "SAMEA13183932", "1000", "75000", "aaa", "100",
        "ftp.example/ERR8972349.fastq.gz")
    paste(paste(.ena_filereport_fields, collapse = "\t"), paste(values, collapse = "\t"), sep = "\n")
}
.v046_mock_ena <- function(perform, sleep = function(seconds) NULL, env = parent.frame()) {
    testthat::local_mocked_bindings(.ena_perform_request = perform,
        .ena_retry_sleep = sleep, .package = "selexprepR", .env = env)
}

test_that("filereport errors with HTTP 200 are retried before parsing", {
    count <- 0L
    delays <- numeric()
    error_text <- paste0(.v046_filereport(), "\nERROR occurred. Not all results may have been written")
    .v046_mock_ena(function(request) {
        count <<- count + 1L
        .v046_ena_response(if (count <= 2L) error_text else .v046_filereport())
    }, function(seconds) delays <<- c(delays, seconds))
    expect_identical(.ena_filereport_text("PRJEB51212", 5), .v046_filereport())
    expect_identical(count, 3L)
    expect_identical(delays, c(1, 2))
})

test_that("persistent filereport failures retain accession and retry diagnostics", {
    count <- 0L
    delays <- numeric()
    .v046_mock_ena(function(request) {
        count <<- count + 1L
        .v046_ena_response("ERROR occurred. Not all results")
    }, function(seconds) delays <<- c(delays, seconds))
    failure <- tryCatch(.ena_filereport_text("PRJEB51212", 5), error = identity)
    expect_s3_class(failure, "error")
    expect_identical(conditionMessage(failure),
        "ENA filereport did not return records after 4 attempts for PRJEB51212; this is the archive's service, not the data: retry later.")
    expect_null(conditionCall(failure))
    expect_identical(count, 4L)
    expect_identical(delays, c(1, 2, 4))
})

test_that("non-success HTTP filereport responses fail before body parsing", {
    count <- 0L
    .v046_mock_ena(function(request) {
        count <<- count + 1L
        .v046_ena_response("ERROR occurred", 400L)
    })
    expect_error(.ena_filereport_text("PRJEB51212", 5), "HTTP status 400")
    expect_identical(count, 1L)
})

test_that("sample XML maps IDs and preserves custom attributes", {
    wanted <- c("SAMEA13183932", "SRS000001", "SAMEA_absent")
    parsed <- .parse_ena_sample_attributes(charToRaw(.v046_sample_xml), wanted)
    expected <- c(construct = "ZFC1 & Sall4", "selex cycle" = "1")
    expect_identical(parsed$SAMEA13183932, expected)
    expect_identical(parsed$SRS000001, expected)
    expect_false("SAMEA_absent" %in% names(parsed))
})

test_that("inspection and download planning use sample rounds and retain evidence", {
    requested <- character()
    timeouts <- numeric()
    .v046_mock_ena(function(request) {
        requested <<- c(requested, request$url)
        timeouts <<- c(timeouts, request$options$timeout_ms)
        if (startsWith(request$url, .ena_filereport_url)) {
            .v046_ena_response(.v046_filereport())
        } else {
            .v046_ena_response(.v046_sample_xml)
        }
    })
    inspection <- inspect_selexprep_accession("PRJEB51212", timeout_seconds = 12)
    result <- fetch_selexprep_reads(inspection, tempfile(), dry_run = TRUE)
    expect_length(requested, 2L)
    expect_identical(timeouts, c(12000, 12000))
    expect_match(requested[[2L]], "SAMEA13183932$")
    expect_identical(inspection$runs$library_strategy, "SELEX")
    expect_identical(inspection$runs$sample_attributes[[1L]],
        c(construct = "ZFC1 & Sall4", "selex cycle" = "1"))
    expect_identical(result$plan$round_number, 1L)
    expect_identical(result$plan$sample_attributes[[1L]], inspection$runs$sample_attributes[[1L]])
    expect_identical(result$plan$library_strategy, "SELEX")
    expect_identical(result$plan$sample_accession, "SAMEA13183932")
})

test_that("a sample XML failure warns and retains text-field fallback", {
    inspection <- .parse_ena_filereport(.v046_filereport(), "PRJEB51212")
    inspection$runs$sample_title <- "Round 2"
    .v046_mock_ena(function(request) stop("network unavailable"))
    observed_warning <- NULL
    enriched <- withCallingHandlers(
        .enrich_ena_sample_attributes(inspection, 5),
        warning = function(condition) {
            observed_warning <<- condition
            invokeRestart("muffleWarning")
        }
    )
    expect_s3_class(observed_warning, "warning")
    expect_identical(conditionMessage(observed_warning),
        "Sample attributes of PRJEB51212 could not be fetched (network unavailable); rounds are parsed from titles and library names only.")
    expect_null(conditionCall(observed_warning))
    expect_identical(enriched$runs$sample_attributes[[1L]], character())
    expect_identical(.fetch_plan(enriched)$round_number, 2L)
})

test_that("malformed sample XML falls back with a warning", {
    inspection <- .parse_ena_filereport(.v046_filereport(), "PRJEB51212")
    .v046_mock_ena(function(request) .v046_ena_response("<SAMPLE_SET><SAMPLE>"))
    expect_warning(enriched <- .enrich_ena_sample_attributes(inspection, 5), "could not be fetched")
    expect_identical(enriched$runs$sample_attributes[[1L]], character())
})

test_that("a missing sample in a batch leaves other sample attributes available", {
    requested <- character()
    .v046_mock_ena(function(request) {
        requested <<- c(requested, request$url)
        if (grepl("SAMD00000001", request$url, fixed = TRUE)) {
            .v046_ena_response(status = 404L)
        } else {
            .v046_ena_response(.v046_sample_xml)
        }
    })
    result <- .query_ena_sample_attributes(c("SAMD00000001", "SAMEA13183932"), 5)
    expect_length(requested, 3L)
    expect_identical(names(result), "SAMEA13183932")
    expect_identical(result$SAMEA13183932[["selex cycle"]], "1")
})

test_that("sample requests are sorted, deduplicated and batched", {
    requested <- character()
    .v046_mock_ena(function(request) {
        requested <<- c(requested, request$url)
        .v046_ena_response("<SAMPLE_SET/>")
    })
    samples <- sprintf("SAMEA%03d", seq_len(51L))
    result <- .query_ena_sample_attributes(c(rev(samples), samples[[1L]], " ", NA_character_), 5)
    expect_length(requested, 2L)
    expect_identical(result, list())
    expect_match(requested[[1L]], "SAMEA001,SAMEA002")
    expect_match(requested[[2L]], "SAMEA051$")
})

test_that("multiple BioSamples assigned to a round are reported without reassignment", {
    inspection <- .parse_ena_filereport(.v046_filereport(), "PRJEB51212")
    second <- inspection$runs
    second$run_accession <- "ERR8972350"
    second$sample_accession <- "SAMEA2"
    second$fastq_urls <- I(list("ftp.example/ERR8972350.fastq.gz"))
    inspection$runs <- rbind(inspection$runs, second)
    inspection$runs$sample_title <- c("Round 4", "Round 4")
    expect_warning(result <- fetch_selexprep_reads(inspection, tempfile()), "more than one BioSample")
    expect_identical(result$plan$round_number, c(4L, 4L))
    expect_identical(S4Vectors::metadata(result$plan)$rounds_with_several_samples,
        list(`4` = c("ERR8972349", "ERR8972350")))
    expect_identical(result$rounds_with_several_samples,
        S4Vectors::metadata(result$plan)$rounds_with_several_samples)

    inspection$runs$sample_accession <- rep("SAMEA1", 2L)
    expect_no_warning(result <- fetch_selexprep_reads(inspection, tempfile()))
    expect_identical(S4Vectors::metadata(result$plan)$rounds_with_several_samples, list())
})

test_that("unknown sample IDs remain conservative and unsafe rounds are excluded from pooling", {
    runs <- S4Vectors::DataFrame(run_accession = c("SRR1", "SRR2", "SRR3"),
        sample_accession = c("", "", "SAMEA3"), sample_title = c("Round 4", "Round 4", "Round 4 / Cycle 7"))
    assignments <- infer_selexprep_rounds(runs)
    expect_identical(.rounds_with_several_samples(runs, assignments), list(`4` = c("SRR1", "SRR2")))
})
