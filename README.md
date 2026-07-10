# selexprep

`selexprep` is an R/Bioconductor package for reproducible preprocessing of
high-throughput SELEX reads. It infers constant regions, extracts variable
regions, creates sparse multi-round count experiments, records provenance, and
reports quality-control signals.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("selexprep")
```

The main workflow is deliberately small:

```r
library(selexprep)

reads_by_round <- list(
    round_00 = read_selexprep_fastq("round_00.fastq.gz"),
    round_01 = read_selexprep_fastq("round_01.fastq.gz")
)
experiment <- run_selexprep(reads_by_round)
```

`experiment` is a `SummarizedExperiment` with a sparse `counts` assay. Its
metadata retains the validated library report, extraction summary, QC report,
and a versioned reproducibility manifest. See the package vignette for the full
workflow.

The bundled study catalog works offline through `selexprep_catalog()`. ENA
requests are explicit: use `inspect_selexprep_accession()` before calling
`fetch_selexprep_reads()`; the latter defaults to a checksum-aware dry-run
plan.

## Development status

This is the `0.99.0` submission candidate for Bioconductor. Development of the
R package is hosted at <https://github.com/marcorotanegroni/selexprepR>; the
original Python implementation remains available at
<https://github.com/marcorotanegroni/selexprep>.
