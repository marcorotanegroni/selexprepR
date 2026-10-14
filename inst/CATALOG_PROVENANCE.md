# Catalog provenance

`selexprep_public_catalog` is a frozen flat view of 240 public HT-SELEX
deposits. It was migrated without modification from
`src/selexprep/catalog/data/curated_metadata.csv` in the original Python
repository at commit `63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5`.

- Snapshot version: `selexprep-v0.4.6`
- Snapshot date: 2026-10-02
- CSV SHA-256:
  `121b599786a718d877bd950660a3a4e0ec178db2f1b54636815d575185c3e0fe`
- Immutable CSV: <https://github.com/marcorotanegroni/selexprep/blob/63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5/src/selexprep/catalog/data/curated_metadata.csv>
- Canonical provenance JSON: <https://github.com/marcorotanegroni/selexprep/blob/63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5/src/selexprep/catalog/data/curated_metadata.json>
- Curation materials: <https://github.com/marcorotanegroni/selexprep/tree/63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5/benchmarks/dual_extraction>

Relative to the v0.4.1 flat catalog, v0.4.6 changes only PRJEB114397
(PFOA)'s target class to `small molecule`; the 240 identifiers and curation
statuses remain unchanged. Its canonical JSON also updates evidence and
concordance notes for PRJEB114397 and PRJNA1481083.

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
