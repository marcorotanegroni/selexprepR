.manifest_synthetic_report <- function() {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    bases <- c("A", "C", "G", "T")
    pools <- stats::setNames(lapply(0:2, function(round) {
        vapply(0:499, function(index) {
            random_region <- paste0(
                bases[((index + round * 1000L) * 7L + 0:19 * 13L) %% 4L + 1L],
                collapse = ""
            )
            paste0(primer_5p, random_region, primer_3p)
        }, character(1))
    }), sprintf("round_%02d", 0:2))
    selexprep_detect(pools)
}

test_that("R-native manifests are deterministic, validated and preserve file hashes", {
    report <- .manifest_synthetic_report()
    input <- tempfile(fileext = ".fastq")
    output <- tempfile(fileext = ".tsv")
    writeLines(c("@read", "ACGT", "+", "IIII"), input)
    writeLines(c("sequence\treads\trank\trpm", "ACGT\t1\t1\t1000000"), output)

    manifest <- build_selexprep_manifest(
        report,
        input_paths = input,
        output_paths = output,
        accession = "PRJNA000001",
        runs = c("SRR000001"),
        parameters = c(min_length = "15"),
        runtime_seconds_per_stage = c(detect = 0.5),
        flags = c("LOW_DEPTH")
    )
    path <- tempfile(fileext = ".json")
    first_hash <- write_selexprep_manifest(manifest, path)
    restored <- read_selexprep_manifest(path)
    second_path <- tempfile(fileext = ".json")
    second_hash <- write_selexprep_manifest(restored, second_path)

    expect_identical(restored$manifest_version, "selexprep_manifest_v2")
    expect_identical(restored$library_report$primer_5p, report$primer_5p)
    expect_identical(restored$input_sha256, manifest$input_sha256)
    expect_identical(restored$output_sha256, manifest$output_sha256)
    expect_identical(readLines(path), readLines(second_path))
    expect_identical(first_hash, second_hash)
    expect_false(endsWith(readChar(path, nchars = file.info(path)$size), "\n\n"))
})

test_that("the manifest reader accepts the portable Python v1 field contract", {
    report <- .manifest_synthetic_report()
    payload <- list(
        manifest_version = "selexprep_manifest_v1",
        selexprep_version = "0.1.0",
        python_version = "3.13.0",
        cutadapt_version = "5.0",
        dnaio_version = "1.2.3",
        pyarrow_version = "18.0.0",
        accession = "PRJNA000001",
        bioproject_id = "PRJNA000001",
        runs = list("SRR000001"),
        input_sha256 = list(reads = "deadbeef"),
        output_sha256 = list(extracted = "cafebabe"),
        library_report = .library_report_payload(report),
        extraction_mode = report$extraction_mode,
        read_source = report$read_source,
        required_action = report$required_action,
        full_insert_recovered = report$full_insert_recovered,
        parameters = list(max_reads = "1000"),
        runtime_seconds_per_stage = list(detect = 1.5),
        flags = list("LOW_DEPTH"),
        sampling_seed = report$sampling_seed
    )
    path <- tempfile(fileext = ".json")
    writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, pretty = TRUE, null = "null"), path)

    manifest <- read_selexprep_manifest(path)

    expect_identical(manifest$manifest_version, "selexprep_manifest_v1")
    expect_identical(manifest$python_version, "3.13.0")
    expect_identical(manifest$library_report$primer_3p, report$primer_3p)
})
