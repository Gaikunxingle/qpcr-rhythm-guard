# Methods and limits

## Normalization

For each sample and gene, arithmetic mean Ct is calculated over observed technical replicates. Missing values are not zero-filled or interpolated. The technical Ct SD is a QC indicator, not a biological standard error. Technical wells do not increase the model sample size.

`delta_ct = mean_ct(target) - mean_ct(reference)`

For each target, the calibrator is the mean complete delta Ct among biological samples in the configured group and elapsed time. It is one constant per target for the whole invocation.

`log2rq = calibrator_mean_delta_ct - delta_ct`

`rq = 2 ** log2rq`

This transformation assumes the relevant equal/doubling amplification-efficiency approximation and a validated stable reference. The software cannot verify these assumptions from Ct alone. It does not propagate calibrator uncertainty into an RQ confidence interval. With an intercept in the regression, subtraction of a common calibrator constant does not change the zero-amplitude test.

## Fixed-period model

The model is ordinary least squares with equal weight per biological sample:

`y = M + C*cos(2*pi*t/P) + S*sin(2*pi*t/P) + error`, with `y = log2rq`.

The period P is fixed in advance (default 24 h). The tool does not estimate or select a period. Amplitude is `sqrt(C*C + S*S)` in log2RQ units. The peak is `(atan2(S,C) * P/(2*pi)) mod P`, relative to input time zero. A numerically zero amplitude has no defined peak. ZT or clock time requires explicit experimental metadata outside this version.

## Inference policy

Descriptive fits can be available when inferential tests are withheld. The ordinary zero-amplitude F test compares an intercept-only model with intercept + cosine + sine, with degrees of freedom 2 and n−3. It tests whether both harmonic coefficients are zero under the fixed period and classical regression assumptions.

The first version permits inference only when all configured and observed requirements pass:

- `design` is explicitly `independent`.
- Each observed sampling time has at least two independent subjects for that gene/group.
- At least four distinct phases modulo P and a time span of at least 0.75 P are available.
- One acquisition run is present; subjects are not reused across time.
- The design matrix is full rank and numerically usable, with positive residual degrees of freedom.
- Target/reference technical SD and calibrator QC permit inference.

These are conservative software policy gates, not universal study-design standards or MIQE recommendations. They do not prove power, independence, normality, homoscedasticity, reference stability, or a biological oscillator. Near-perfect fits without meaningful residual variance are not assigned spurious tiny P values. Constant responses have amplitude zero and, when otherwise eligible, P=1. R-squared is undefined when the response is constant.

Missing normalized rows remain in the output. Available-case fitting may be biased when nondetection depends on expression. The user must inspect missingness and residual behavior before interpreting a test. A threshold warning retains descriptive values; it does not silently remove inconvenient observations.

## Multiple tests and unsupported claims

Benjamini–Hochberg adjustment is applied to all finite P values for all gene-by-group models in this invocation. Untestable models retain NA and are not treated as P=1. The manifest records planned/tested model counts. This family does not automatically include tests performed in other runs; splitting runs can change the family and requires a predeclared analysis plan.

There are no between-group waveform, amplitude or phase tests; no repeated-measure or batch-adjusted inference; no endogenous rhythm certification. “Significant in one group and nonsignificant in another” is not evidence of a group difference.

## Sources and related software

- Cornelissen G. *Cosinor-based rhythmometry.* 2014. [Full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3991883/). Model structure, study designs and regression assumptions; the software policy thresholds above are our implementation choices.
- Bustin et al. *MIQE 2.0: Revision of the Minimum Information for Publication of Quantitative Real-Time PCR Experiments Guidelines.* 2025. [PubMed record](https://pubmed.ncbi.nlm.nih.gov/40272429/). Reporting context; this project does not implement or certify the full guideline.
- *CosinorPy: a python package for cosinor-based rhythmometry.* [Article DOI](https://doi.org/10.1186/s12859-020-03830-w). A related established package. Our prototype focuses on the Ct-to-QC audit trail; it does not claim to introduce cosinor analysis or replace broader rhythmometry software.
