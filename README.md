# selexprepR

`selexprepR` is an R/Bioconductor package for reproducible preprocessing of
high-throughput SELEX reads. It infers constant regions, extracts variable
regions, creates sparse multi-round count experiments, records provenance, and
reports quality-control signals.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("selexprepR")
```

The main workflow is deliberately small:

```r
library(selexprepR)

experiment <- run_selexprep_files(c(
    round_00 = "round_00.fastq.gz",
    round_01 = "round_01.fastq.gz"
))
```

Inline-barcoded pools can be split with an explicit barcode-to-round map:

```r
inputs <- selexprep_demultiplex(
    read_selexprep_fastq("multiplexed.fastq.gz"),
    c(AAAAA = 0L, TTTTT = 1L)
)
experiment <- run_selexprep(inputs)
```

`experiment` is a `SummarizedExperiment` with a sparse `counts` assay. Its
metadata retains the validated library report, extraction summary, QC report,
and a versioned reproducibility manifest. See the package vignette for the full
workflow.

The bundled study catalog is both a package dataset and a queryable table:

```r
data("selexprep_public_catalog", package = "selexprepR")
selexprep_catalog("FGF-9")
```

It works offline and carries snapshot provenance in
`S4Vectors::metadata(selexprep_public_catalog)`. ENA requests are explicit: use
`inspect_selexprep_accession()` before calling `fetch_selexprep_reads()`; the
latter defaults to a checksum-aware dry-run plan.

The R implementation is aligned to Python `selexprep v0.4.6` (commit
`63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5`), with independent R versioning.
Shared constants are checked across rounds; extraction warns when fewer than
half the reads survive, with separate FASTQ diagnostics for file workflows. `HIGH` and `MEDIUM` describe read support rather than a
verified constant/random boundary. Consult the study's library design when a
diverse early round is missing. See `inst/UPSTREAM_ALIGNMENT.md` for scope,
intentional R interface differences, and reproducible regression references.

## Development status

This is the `0.99.6` submission candidate for Bioconductor. Development of the
R package is hosted at <https://github.com/marcorotanegroni/selexprepR>; the
original Python implementation remains available at
<https://github.com/marcorotanegroni/selexprep>.
