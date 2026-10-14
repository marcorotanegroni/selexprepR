# selexprepR 0.99.6

* Simplify ENA condition messages to address the two new BiocCheck coding
  practice notes, preserving four-attempt retry and sample-metadata fallback.
* Keep diagnostic accessions, underlying sample failures and retry advice
  explicit; remove redundant failure wording from the filereport message.

# selexprepR 0.99.5

* Align library inference with Python selexprep v0.4.6: retain shared library
  constants across rounds using per-round offsets, preserve the insert boundary,
  and avoid ambiguous repeated-core anchors.
* Warn about extraction pools retaining fewer than half their reads or half
  the best pool's yield, retaining diagnostics in `low_yield_rounds`. File
  workflows also retain `input_yields` and `low_yield_inputs` per FASTQ, so
  aggregation cannot mask a poorly matching input.
* Normalize lowercase RNA consistently during detection and primer matching.
* Read ENA BioSample attributes for round inference; handle conflicting and
  weak PCR-cycle evidence, and stop interpreting ArrayExpress RV sample labels
  as selection rounds.
* Retry ENA filereport service-error bodies, retain run-level library strategy
  and sample evidence, and warn when a round pools distinct BioSamples.
* Pin the 240-record catalog and provenance to Python v0.4.6, including the
  corrected PFOA target class. Add offline Python/R regression references.

# selexprepR 0.99.4

* Refuse inferred libraries with a zero-length random region instead of
  reporting a valid full-insert extraction.
* Rescue well-supported primer cores when noisy outer flanks lower the
  full-length match rate, without changing the extraction boundary.
* Update the bundled public catalog to the adjudicated Python v0.4.1
  snapshot, retaining immutable source and evidence provenance.

# selexprepR 0.99.3

* Write LibraryReport and manifest JSON with UTF-8 LF line endings on every
  platform, preserving deterministic hashes and Python golden-file parity on
  Windows.

# selexprepR 0.99.1

* Initial Bioconductor submission candidate.
* Provides primer inference, extraction, sparse multi-round counting, quality
  control, reproducibility manifests, FASTQ input, and ENA data discovery.
* Adds `selexprep_public_catalog` as a documented, reproducibly generated
  package dataset with immutable snapshot provenance.
* Adds R-native barcode demultiplexing and direct FASTQ/fetch-result input
  bridges with file hashes retained in the run manifest.
