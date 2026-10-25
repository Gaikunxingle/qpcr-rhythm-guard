# Contributing

Start with a small reproducible example and the expected behavior. Include the software version, configuration, relevant QC output, and whether the data are synthetic. Do not post unpublished patient/sample records, credentials or restricted laboratory files.

Run `python -m unittest discover -s tests -v` and the synthetic CLI example before proposing a change. New analytical behavior needs an independent numerical reference and an explicit description of its assumptions. Preserve missing/untestable outputs and label provenance.

For bug reports, state the input, observed output, expected output and environment. For features, explain the real workflow that is blocked. The maintainer is [Gaikunxingle](https://github.com/Gaikunxingle); use [repository Issues](https://github.com/Gaikunxingle/qpcr-rhythm-guard/issues) for public questions and reproducible reports.
