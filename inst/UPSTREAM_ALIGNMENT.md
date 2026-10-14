# Alignment to Python selexprep v0.4.6

R 0.99.5 aligns with the Python release tag v0.4.6, commit
`63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5` (2026-10-02). The prior verified
reference was v0.4.1, commit `fdb110736365261c5f000050dc77f08d41fa33fd`.
Code, tests, catalog CSV/JSON and documentation were compared between those
commits. Later Python main development is outside this alignment.

## Runtime changes

| Area | Behavior aligned |
| --- | --- |
| Library inference | Shared constants retain the inner 12 nt at the insert boundary. Only outer bases unsupported by anchored rounds are removed; each round uses its own offset for positional signals. Single-copy agreement resolves duplicate cores, with conservative fallback for ambiguous flanks. Full earliest-round flanks remain the basis of insert length and orientation. |
| Extraction | Nonempty pools below 50% own retention or half the best yield warn and retain evidence. File workflows measure each R1 FASTQ separately, with joint retention of split-primer pairs, before aggregation can hide a failed input. Empty files have an undefined fraction and do not warn. |
| Alphabet | Detection and primer matching uppercase RNA before converting U to T, including lowercase u. Extracted insert case/alphabet remain unchanged as in the R API. |
| ENA | HTTP-200 filereport bodies containing `ERROR occurred` retry four times with 1/2/4-second delays. Custom sample XML attributes are fetched in sorted, deduplicated batches of 50; missing individual records are omitted, and sample-service failure warns and falls back to text. Run-level strategy and sample evidence are retained. |
| Rounds | RV sample labels no longer assign rounds. Normalized structured keys and prefixed integer values are recognized; strong conflicting attributes remain unsafe. PCR-cycle evidence is MEDIUM and does not override explicit text evidence. Distinct BioSamples sharing an assigned round warn without changing assignments. |
| Catalog | The same 240 identifiers and schema are retained; PRJEB114397 (PFOA) is classified as `small molecule`. CSV and evidence provenance are pinned to the release commit; see CATALOG_PROVENANCE.md. |

The pre-existing zero-insert guard, noisy-flank core rescue and LibraryReport
JSON schema already matched v0.4.1 and are preserved. HIGH/MEDIUM measure read
support, not a verified constant/random boundary; consult the library design
when diverse early rounds are absent.

## R interface choices retained

- Named round pools, S3 decision records, DataFrames and a sparse
  SummarizedExperiment remain the R interface. Missing values use NA/NULL.
- R extraction uses Biostrings rather than a cutadapt subprocess. Agreement is
  tested on shared fixtures; arbitrary approximate-match tie cases are not
  claimed to be exhaustively equivalent.
- ENA inspection is explicit; fetch defaults to dry-run and refuses the whole
  unsafe plan. Manual overrides remain named vectors. No external download
  backend or CLI-only API is added.
- R records both candidates in `Round 1 Round 2`. Python v0.4.6 records only
  the first because its repeated regex consumes the shared separator. R keeps
  its conservative existing ambiguity handling.
- Round/FASTQ diagnostics are named R vectors and DataFrames; CLI logging,
  exact message formatting and benchmark tables are not public R contracts.

Python benchmark metrics, simulation tooling and generated evaluation results
are outside this package alignment. Counting, QC and manifest schemas do not
gain unrelated features.

## Reproducible regression evidence

All new package tests run offline and require no Python installation:

- `v046-detection-reference.json`: all stable fields of eight complete
  LibraryReports from identical 600-read pools, including variable tags,
  repeated insert/flank cores, heterogeneous outer bases and paired reads.
  Compact pool construction and exact source provenance are documented beside
  the fixture in `v046-detection-reference-provenance.md`.
- `v046-extraction-reference.json`: 108 identical reads across three cases
  and three rounds, extracted by the tagged Python trim wrapper with cutadapt;
  R compares every emitted insert. Inputs are the first 12 reads of each
  corresponding detection-reference round, using its full reference report.
- `v046-round-reference.json`: 41 metadata inputs executed with
  `parse_round` extracted directly from the tag, compared for round number,
  confidence, evidence source, candidates and assignment safety.
- Focused regressions cover exact support/coverage/anchor thresholds,
  overlapping cores, service-error retry, XML aliases/batching/fallback,
  attribute precedence, distinct BioSamples, low-yield file aggregation,
  paired retention and empty/capped files.

An additional independent audit compared 350 generated shared-constant helper
cases to functions extracted from the exact Python tag; all agreed. This is
supplemental audit evidence, not a replacement for the committed offline tests.

R version numbering is independent of Python. The requirement R >= 4.6.0
is unchanged. Local source-mode tests on older R do not substitute for the
configured R 4.6.1/Bioconductor devel build, vignette, check and BiocCheck CI.
