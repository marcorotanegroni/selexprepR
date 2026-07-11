# AI-assisted development disclosure

The migration from the original Python implementation to this R package used
OpenAI Codex as a development assistant for code translation, test design,
documentation, and Bioconductor compliance review. The bundled catalog has a
separate, more specific dual-model curation history documented in
`CATALOG_PROVENANCE.md`.

The maintainer reviewed the implementation against the original package
contracts and remains responsible for correctness, licensing, security,
maintenance, and user support. Stable JSON fixtures produced by the Python
implementation, focused R unit tests, and differential extraction cases against
cutadapt are used to detect behavioral drift. No Python runtime or generated
foreign-language source code is bundled in the R package.

This disclosure should also be included in the BiocContributions submission
issue, together with any additional AI assistance used after this snapshot.
