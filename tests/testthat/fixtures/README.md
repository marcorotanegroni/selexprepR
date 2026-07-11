# Python parity fixtures

These `LibraryReport` JSON files were produced from the original Python
implementation at commit `b6792637ec9bed78f131e8f63e12e155f733e9ae` using the
deterministic synthetic pools reconstructed in `test-library-report.R`. They
cover both-primer, 5-prime-only, 3-prime-only, paired split-primer, and
safe-failure classifications.

The R tests validate both parsed field equality and byte-for-byte JSON output.
Reviewed fixture SHA-256 values are:

```text
e6cc473e74fb0269eb779b93871bfdccff0ee943fedd7ef7ab49ef00597ac627  library-report-both-primers.json
47035ac609795c2acbc77cab1629644780d27c70f5d21eb9597f050ccc7062cd  library-report-five-prime-only.json
d5d499ce4f1c39bf2310db39b954804e0028155d50d4b3439238abef73acd17b  library-report-paired-split.json
0c817aa1ae705a9e2c4af7dfca3575ab893ea317ab7138970f36a89e29b298fc  library-report-three-prime-only.json
5cf6c206cc65bf3ba212a86a53e36ec2e9cb0ffaa4da83ac1a77eb02d0d726fe  library-report-unable.json
```

Extraction parity is tested separately against cutadapt 5.2 for substitutions,
insertions, deletions, unanchored full-primer matches, RNA/DNA alphabet
normalization, and paired-end coordinate handling.
