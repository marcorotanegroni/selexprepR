# Catalog provenance

`selexprep_public_catalog` is a frozen flat view of 240 public HT-SELEX
deposits. It was migrated without modification from
`src/selexprep/catalog/data/curated_metadata.csv` in the original Python
repository at commit `b6792637ec9bed78f131e8f63e12e155f733e9ae`.

- Snapshot version: `v0.2.1-dual-extraction-2026-07-03`
- Snapshot date: 2026-07-03
- CSV SHA-256:
  `9e11bc3868816e769da3ae50c159bb14e679e790be6de106b0ede3756caa1cb9`
- Immutable CSV: <https://github.com/marcorotanegroni/selexprep/blob/b6792637ec9bed78f131e8f63e12e155f733e9ae/src/selexprep/catalog/data/curated_metadata.csv>
- Canonical provenance JSON: <https://github.com/marcorotanegroni/selexprep/blob/b6792637ec9bed78f131e8f63e12e155f733e9ae/src/selexprep/catalog/data/curated_metadata.json>
- Curation materials: <https://github.com/marcorotanegroni/selexprep/tree/b6792637ec9bed78f131e8f63e12e155f733e9ae/benchmarks/dual_extraction>

The eight experimental fields were extracted independently by Claude and
Codex/GPT, then reconciled. The paired `_curation` columns preserve whether a
cell was concordant, discordant, absent, independently verified, or supplied by
one extraction arm. Discordant values retain both outputs, separated by
`" || "`. Evidence quotations are not redistributed in the R package; the
immutable canonical JSON above preserves record-level sources and locations.

The package code is distributed under the MIT license. Catalog values are
factual metadata assembled from public repository records and cited studies;
users should consult the terms of the originating record or publication before
redistributing source material. Inclusion in this snapshot does not assert
that the corresponding reads remain available.
