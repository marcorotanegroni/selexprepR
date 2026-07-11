test_that("selexprep_qc computes per-round diversity summaries", {
    experiment <- selexprepR:::.as_selexprep_experiment(list(
        round_00 = selexprep_count(c("AAAA", "AAAA", "CCCC", "GGGG")),
        round_01 = selexprep_count(c("AAAA", "AAAA", "AAAA", "AAAA"))
    ))
    qc <- selexprep_qc(experiment, low_total_reads = 0)

    expect_s3_class(qc, "selexprep_qc")
    expect_true(all(c(
        "round", "n_reads", "n_unique", "shannon_entropy_bits",
        "rarefied_unique", "singleton_fraction", "top_1_coverage",
        "top_100_coverage", "modal_sequence_length", "mean_gc", "gc_sd",
        "nonstandard_alphabet_fraction", "truseq_reads_fraction", "top_kmer_fraction"
    ) %in% colnames(qc$per_round)))
    expect_equal(qc$per_round$n_reads, c(4, 4))
    expect_equal(qc$per_round$n_unique, c(3L, 1L))
    expect_equal(qc$per_round$shannon_entropy_bits, c(1.5, 0))
    expect_equal(qc$per_round$top_100_coverage, c(1, 1))
})

test_that("selexprep_qc reports low depth and report-derived flags", {
    experiment <- selexprepR:::.as_selexprep_experiment(list(
        round_00 = selexprep_count(c("AAAA", "CCCC"))
    ))
    fields <- unclass(read_library_report(test_path("fixtures", "library-report-both-primers.json")))
    fields$match_rate_5p <- 0.1
    fields$required_action <- "READ_MERGING_RECOMMENDED"
    fields$extraction_mode <- "PAIRED_END_SPLIT_PRIMERS"
    report <- selexprepR:::.new_library_report(fields)

    qc <- selexprep_qc(experiment, report, low_total_reads = 10)
    expect_setequal(qc$flags$name, c(
        "low_total_reads", "low_primer_match", "requires_read_merging_for_full_insert"
    ))
})

test_that("QC rarefaction is deterministic and surfaces diversity increases", {
    counts <- c(A = 100, B = 50, C = 25)
    sampled_a <- selexprepR:::.rarefy_counts(counts, depth = 50, seed = 42)
    sampled_b <- selexprepR:::.rarefy_counts(counts, depth = 50, seed = 42)
    expect_identical(sampled_a, sampled_b)
    expect_equal(sum(sampled_a), 50)
    expect_true(all(sampled_a > 0))

    bases <- c("A", "C", "G", "T")
    diverse <- vapply(0:9, function(index) {
        paste0(bases[(index %/% 4^(0:5)) %% 4 + 1], collapse = "")
    }, character(1))
    experiment <- selexprepR:::.as_selexprep_experiment(list(
        round_00 = selexprep_count(rep("AAAAAA", 10)),
        round_01 = selexprep_count(diverse)
    ))
    qc <- selexprep_qc(experiment, low_total_reads = 0, rarefaction_depth = 10)

    expect_setequal(qc$flags$name, "unexpected_rarefied_diversity_increase")
    expect_equal(qc$per_round$rarefied_unique, c(1L, 10L))
    expect_identical(qc$settings$effective_rarefaction_depth, 10)
})

test_that("QC identifies sequence-quality and adapter diagnostics", {
    experiment <- selexprepR:::.as_selexprep_experiment(list(
        round_00 = selexprep_count(c(
            rep("AGATCGGAAGAGC", 20), "ACGTAC", "ACGTXG"
        ))
    ))
    fields <- unclass(read_library_report(test_path("fixtures", "library-report-both-primers.json")))
    fields$known_adapter_hits <- c(TRUSEQ_R1 = 10, NEXTERA = 0)
    report <- selexprepR:::.new_library_report(fields)
    qc <- selexprep_qc(experiment, report, low_total_reads = 0)

    expect_true(all(c("nonstandard_alphabet", "truseq_residual", "adapter_contamination_high") %in%
        qc$flags$name))
    adapter_row <- which(qc$flags$name == "adapter_contamination_high")
    expect_identical(qc$flags$severity[[adapter_row]], "info")
})

test_that("QC flags k-mer progression only as a diagnostic", {
    distances <- S4Vectors::DataFrame(
        from_round = c("round_00", "round_00", "round_01"),
        to_round = c("round_01", "round_02", "round_02"),
        jaccard_distance = c(0.8, 0.1, 0.2)
    )
    flag <- selexprepR:::.round_kmer_monotonicity_flag(distances)

    expect_identical(flag$name, "round_kmer_nonmonotonic")
    expect_identical(flag$severity, "warn")
    expect_length(flag$evidence$violations, 1L)
})

test_that("QC surfaces a high reverse-strand fraction from extraction provenance", {
    primer_5p <- "GGTAATACGACTCACTATAGGG"
    primer_3p <- "CCATGCATGCATGCATGCAT"
    fields <- unclass(read_library_report(test_path("fixtures", "library-report-both-primers.json")))
    fields$orientation <- "MIXED"
    report <- selexprepR:::.new_library_report(fields)
    forward <- paste0(primer_5p, "ACGT", primer_3p)
    extraction <- selexprep_extract(
        list(round_00 = c(forward, reverse_complement(forward))),
        report
    )
    experiment <- selexprepR:::.as_selexprep_experiment(list(
        round_00 = selexprep_count(extraction$sequences_by_round$round_00)
    ))
    metadata <- S4Vectors::metadata(experiment)
    metadata$extraction <- extraction
    metadata$library_report <- report
    S4Vectors::metadata(experiment) <- metadata

    qc <- selexprep_qc(experiment, low_total_reads = 0)

    expect_true("strand_mix" %in% qc$flags$name)
})
