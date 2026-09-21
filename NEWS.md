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
