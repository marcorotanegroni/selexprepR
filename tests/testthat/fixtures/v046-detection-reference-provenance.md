# Python v0.4.6 detection reference

`v046-detection-reference.json` contains complete LibraryReports and the exact
input sequences for eight offline Python/R regression comparisons. The Python
reference is the release tag `v0.4.6`, commit
`63f33b44bb381c2c6b6ed131f3b98e8d6ff9d5f5`, in
<https://github.com/marcorotanegroni/selexprep>. The earlier alignment reference
was `v0.4.1`, commit `fdb110736365261c5f000050dc77f08d41fa33fd`.

Reports were generated with `compute_library_report()` from an isolated copy
of the tagged Python source, using Python 3.12.13, NumPy 2.4.5 and pandas 3.0.3.
Detection used its defaults, including sampling seed 42 and no read cap. No
code or data from subsequent Python development was used.

To avoid repeating reads across scenarios, the fixture stores seven pools of
600 random 40-nt inserts. Pool names are the seeds of Python's
`random.Random(seed)`. Each insert uses 40 consecutive
`rng.choice("ACGT")` calls. Round descriptors specify a literal prefix and
suffix to concatenate with every insert. `leading_cycle: "GAC"` prepends
successive G, A and C bases, repeating from G. `insert_edge_copies: true`
replaces the first and last 12 nt of the insert with `core_5p` and `core_3p`
for zero-based read indices whose remainder modulo 5 is 0, 1 or 2. The
construction was checked against every original uncompressed reference read.
The test helper reconstructs these reads without Python, network access or a
random-number generator.

The constants are `AGGACGATGCGGTACC` (5-prime) and
`CAGACGACTCGCTGAGGATC` (3-prime). The cases cover:

- A clean library across three rounds.
- Different run-specific tags, including unequal tag lengths.
- Those tags with additional core copies inside later-round inserts.
- A 3-prime outer base that differs between rounds.
- A heterogeneous leading base in a single round.
- Repeated core copies in the 5-prime flank.
- Repeated core copies in the 3-prime flank.
- Paired split-primer reads with run-specific tags on R1 and R2.

All stable LibraryReport fields are compared. These fixtures preserve the
distinction between a shared library constant and the complete earliest-round
flank used to determine the original insert boundary and orientation.
