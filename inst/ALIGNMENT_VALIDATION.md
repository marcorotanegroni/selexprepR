# Local validation and submission update: selexprepR 0.99.5

Validated on 2026-10-14 using the repositories at
`/Users/marcorotanegroni/selexprepR` and `/Users/marcorotanegroni/selexprep`.
The R working tree initially had no modifications; Python had only an
untracked `AGENTS.md`, which remains unchanged. R started at 0.99.4,
`a6efea0d92b2e496b334a725396cf2ffedeae4e0`, on main.

The reference is Python v0.4.6,
`63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5`. The comparison baseline is
v0.4.1, `fdb110736365261c5f000050dc77f08d41fa33fd`. Later Python main commits
were excluded. See `UPSTREAM_ALIGNMENT.md` for the runtime changes,
intentional R differences and regression-fixture provenance.

## Results

| Verification | Result | Qualification |
| --- | --- | --- |
| Full R source-mode suite | 95 test blocks; 663 assertions passed; zero failed, skipped, errors or test warnings | R 4.3.1 via testthat/pkgload; this is below the supported R requirement |
| Same-input Python/R LibraryReports | Eight complete reports agree within numerical tolerance | Variable tags, duplicate cores, outer-base differences and paired split |
| Same-input Python/R extraction | Every emitted insert agrees for 108 reads in nine pools | Tagged Python trim wrapper with cutadapt; does not exhaustively establish all approximate-match tie cases |
| Same-input Python/R round inference | 41 metadata cases agree on number, confidence, source, candidates and safety | Executed parse_round source extracted from the exact tag |
| Independent shared-constant review | 350 further generated helper cases agree | Supplemental audit, in addition to committed tests |
| Selected exact-tag Python regressions | 105 tests passed | Six tag-pinned modules; includes shared constants, extraction yield, metadata, ENA retry, fetch plan and status note; Python repository was not modified |
| Catalog | 240 identifiers; expected CSV SHA-256; rebuilt R dataset and metadata match | CSV SHA-256 `121b599786a718d877bd950660a3a4e0ec178db2f1b54636815d575185c3e0fe` |
| Generated documentation | All Rd files parse and check; public exports preserved | Roxygen regenerated relevant help files |
| Source build without vignette rendering | Succeeded, selexprepR_0.99.5.tar.gz | R CMD build --no-build-vignettes; not a full supported-environment build |
| Full source build | Blocked during installation | Local R is 4.3.1; DESCRIPTION requires R >= 4.6.0 |
| R CMD check --no-manual --no-build-vignettes | Stops at dependency check, 1 ERROR | Suggested packages arrow and BiocStyle unavailable; no complete check result |
| BiocCheck and BiocCheckGitClone | Not run locally | BiocCheck is not installed |
| Git whitespace and patch integrity | Passed | No incomplete patches or Python modifications |

The R source tests use:

```r
Sys.unsetenv("LC_ALL")
testthat::test_local("/Users/marcorotanegroni/selexprepR",
                     reporter = "summary", stop_on_failure = TRUE)
```

The installed mixed Bioconductor environment emits a pre-existing namespace
import replacement warning for S4Arrays/DelayedArray while loading
SummarizedExperiment. It is outside the test assertions; the suite itself has
zero warnings. Local Bioconductor is 3.17. BiocStyle, arrow, BiocCheck,
pandoc and pdflatex are unavailable. A check retry with network access removed
the initial sandbox index warnings but still correctly stopped at missing
Suggests. No dependency requirement or check setting was relaxed to obtain a
passing result.

The package retains `Depends: R (>= 4.6.0)` and the author name Marco Rota
Negroni. The configured CI uses R 4.6.1/Bioconductor devel, installs the new
xml2 import explicitly, and runs BiocCheckGitClone, a full build with
vignettes, R CMD check and BiocCheck new-package. Those results remain pending
for this candidate, including the destination submission's platform builds.

## Recovery after the usage-limit interruption

The three subagents had written only temporary copies and patch artifacts.
After resuming, each agent verified its saved state. Completed patches were
reviewed and applied by the primary agent; the per-file diagnostics were
finished before integration. Focused tests were repeated against the actual
R repository, followed by the full suite and package checks above. No partial
subagent edit remained in either repository. The Python working tree still
contains only its original untracked AGENTS.md.

## Existing submission: maintainer's next steps

On 2026-10-14, read-only remote inspection found the personal main and
BiocStaging devel branches both at 0.99.4 and commit a6efea0, without tags
reserving 0.99.5. The existing
[submission #90](https://github.com/Bioconductor/BiocContributions/issues/90)
remains open. Its
[bot instruction](https://github.com/Bioconductor/BiocContributions/issues/90#issuecomment-5775033179)
requires a version increment to produce a new report.

1. Review the local candidate and NEWS, then recheck remote heads for changes
   since this audit. Preserve any new remote work; do not force-push.
2. Publish the prepared main commit to marcorotanegroni/selexprepR and check
   its R 4.6.1/Bioconductor devel CI.
3. Update the existing BiocStaging/selexprepR devel branch with the same
   history, using the existing remote mapping main:devel.
4. Wait for the 0.99.5 report in issue #90; review source, Bioc-checks, Linux,
   macOS, Windows and wasm. Local source-mode tests do not replace these.
5. Continue the existing review and update the AI-assistance disclosure when
   communicating with the reviewer. No new submission or GitHub release is
   needed merely to trigger this report.

No push, release creation, issue comment or remote write was performed.
