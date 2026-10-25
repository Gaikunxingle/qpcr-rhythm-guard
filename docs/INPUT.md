# Input contract

UTF-8 (optional BOM), comma-separated CSV, one row per technical qPCR measurement. All eight named columns are required. Labels are case-sensitive. Blank identifiers are invalid.

| Column | Meaning |
|---|---|
| sample_id | Unique physical sample. Its subject, group, time and run must agree on every row. |
| subject_id | Independent biological unit, e.g. independently randomized animal/culture. |
| group | Exact experimental group label. |
| time_h | Finite, nonnegative elapsed hours from a common declared origin. |
| gene | Assay/target name, including the configured reference gene. |
| tech_rep | Technical replicate identifier unique within sample and gene. |
| ct | Finite Ct number or a documented missing token. Zero is not a missing token. |
| run_id | Acquisition run/batch. Multiple runs disable inference in this version. |

Missing tokens are blank, NA, N/A, Undetermined, Unknown, None and Null (case insensitive). Numeric NaN and Infinity are rejected. Undetected assays may be non-random missingness; missing rows are retained, not imputed.

Every sample is paired with its reference gene within that sample only. The same subject at one time cannot be assigned to multiple sample IDs. A subject cannot belong to multiple groups. Repeated subjects across different times may be represented, but inference is disabled. This narrow format does not cover crossover or multiple-tissue designs.

Configuration keys: `reference_gene`, `calibrator_group`, `calibrator_time_h` (required), `period_h` (default 24), `design` (unknown, independent, repeated; default unknown), `technical_sd_limit` (default 0.5). Typos/unknown keys are rejected. JSON examples are in `examples/`.

For spreadsheet safety, exported CSV text beginning with a formula prefix is escaped with an apostrophe. Exact identifiers are retained in JSON. CSV numeric missing values are empty; the HTML report renders them as NA.
