# qPCR Rhythm Guard

[中文使用说明](README.zh-CN.md)

[Source](https://github.com/Gaikunxingle/qpcr-rhythm-guard) · [Report an issue](https://github.com/Gaikunxingle/qpcr-rhythm-guard/issues)

Local qPCR quality checks, explicit biological-replicate handling, and conservative fixed-period cosinor analysis. An early research prototype, version 0.1.0. No cloud upload, model subscription, or API key is needed to run it.

The practical question is: **which Ct records can be normalized and fitted, and which assumptions or missing observations need attention?** This project combines traceable preprocessing with limited, clearly labelled inference. It is not a general circadian discovery platform or a substitute for assay validation.

## Quick start

Requires Python 3.10 or later. In this project folder:

```sh
python -m pip install .
python -m qpcr_rhythm_guard examples/synthetic.csv --config examples/config.json --out results/demo
```

Open `results/demo/report.html`. The included data are entirely synthetic. `design: unknown` is the safe default: curves are descriptive and P values are not computed. A separate artificial example configuration exercises the independent-sample inference branch:

```sh
python -m qpcr_rhythm_guard examples/synthetic.csv --config examples/config_independent_synthetic.json --out results/synthetic_inference
python -m unittest discover -s tests -v
```

Output folders must be empty or new. Existing results are preserved.

## Implemented

- Read a documented UTF-8 CSV format; reject duplicate records and conflicting sample metadata.
- Retain missing/Undetermined Ct values and reasons; never replace them with zero.
- Average technical Ct replicates within a sample and gene; report technical scatter.
- Normalize with one declared reference gene and a gene-specific group/time calibrator.
- Fit fixed-period cosinor in log2RQ; do not estimate a free-running period.
- Gate the zero-amplitude F test on explicit independent-sample design and conservative input checks.
- Apply BH adjustment to all testable gene-by-group models in one invocation.
- Export HTML, CSV, JSON, PNG/SVG, input/config hashes, and software versions.

## Input

See [the CSV schema](docs/INPUT.md). Required columns:

```text
sample_id,subject_id,group,time_h,gene,tech_rep,ct,run_id
```

`subject_id` identifies the independently randomized biological unit, not a qPCR well. Do not relabel repeated measurements as different subjects to enable a test. Elapsed hours are numeric; sample order and injury stage are not clock time.

## Scientific scope

Read [METHODS](docs/METHODS.md) before interpreting output. Single-reference stability and amplification efficiency are not established from these input Ct values. Raw amplification curves, NTC/no-RT controls, multiple-reference normalization, plate correction, mixed-effects repeated-measure models, between-group rhythm tests, native Prism export, and instrument-specific Excel/EDS import are not implemented.

The application thresholds are software guardrails, not MIQE certification or proof of adequate power. P values do not establish an endogenous circadian mechanism. Missing or untestable outcomes remain in the output.

## Development and maintenance

See [CONTRIBUTING](CONTRIBUTING.md) and [the roadmap](docs/ROADMAP.md). [GitHub Actions](https://github.com/Gaikunxingle/qpcr-rhythm-guard/actions/workflows/tests.yml) checks installation, 22 tests, and synthetic report generation on Windows and Ubuntu with Python 3.10 and 3.12. All four jobs passed in the [initial publication run](https://github.com/Gaikunxingle/qpcr-rhythm-guard/actions/runs/40176021104) on 2026-10-25. No external adoption, publication, download count, or grant award is claimed.

MIT licensed. Maintained by [Gaikunxingle](https://github.com/Gaikunxingle). The package-index name has not been reserved; install from this repository.
