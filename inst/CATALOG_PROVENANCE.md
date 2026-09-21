# Catalog provenance

`selexprep_public_catalog` is a frozen flat view of 240 public HT-SELEX
deposits. It was migrated without modification from
`src/selexprep/catalog/data/curated_metadata.csv` in the original Python
repository at commit `0db8e3313a775e30c7bc961aa61d121e75c74467`.

- Snapshot version: `v0.3.2-dual-extraction-adjudicated-en-2026-09-04`
- Snapshot date: 2026-09-04
- CSV SHA-256:
  `e67b7b77f9d60b6ae6a686099c2a524ae59f50fbaa446aece7ea607fdc007fef`
- Immutable CSV: <https://github.com/marcorotanegroni/selexprep/blob/0db8e3313a775e30c7bc961aa61d121e75c74467/src/selexprep/catalog/data/curated_metadata.csv>
- Canonical provenance JSON: <https://github.com/marcorotanegroni/selexprep/blob/0db8e3313a775e30c7bc961aa61d121e75c74467/src/selexprep/catalog/data/curated_metadata.json>
- Curation materials: <https://github.com/marcorotanegroni/selexprep/tree/0db8e3313a775e30c7bc961aa61d121e75c74467/benchmarks/dual_extraction>

The eight experimental fields were extracted independently by Claude and
Codex/GPT, then reconciled. The paired `_curation` columns preserve whether a
cell was concordant, adjudicated after disagreement, absent, independently
verified, or supplied by one extraction arm. The 47 adjudicated cells have a
single resolved flat value; the canonical JSON preserves both extraction arms,
the resolution rationale, and record-level evidence sources and locations.

The package code is distributed under the MIT license. Catalog values are
factual metadata assembled from public repository records and cited studies;
users should consult the terms of the originating record or publication before
redistributing source material. Inclusion in this snapshot does not assert
that the corresponding reads remain available.
