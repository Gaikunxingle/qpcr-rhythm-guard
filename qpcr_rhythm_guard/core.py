"""qPCR normalization and fixed-period cosinor with conservative testing gates.

``analyze`` preserves input identifiers and raw CSV values. A sample is the
biological observation; technical wells are aggregated before normalization.
The single reference gene and equal amplification-efficiency assumptions are
specified by the caller and cannot be validated from Ct values alone.
"""

from __future__ import annotations

import csv
import hashlib
import io
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import f as f_distribution


REQUIRED_COLUMNS = (
    "sample_id", "subject_id", "group", "time_h", "gene", "tech_rep", "ct", "run_id"
)
MISSING_CT_TOKENS = frozenset({"", "na", "n/a", "undetermined", "unknown", "none", "null"})
CONFIG_KEYS = frozenset({
    "reference_gene", "calibrator_group", "calibrator_time_h", "period_h",
    "design", "technical_sd_limit",
})


def _finite_number(value: Any, name: str, *, minimum: float | None = None,
                   positive: bool = False) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number, not a boolean")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite; NaN and infinity are not accepted")
    if minimum is not None and result < minimum:
        raise ValueError(f"{name} must be >= {minimum:g}")
    if positive and result <= 0:
        raise ValueError(f"{name} must be > 0")
    return result


def _identifier(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    # Deliberately preserve leading/trailing spaces and case in identifiers.
    return value


def _validate_config(config: dict[str, Any]) -> dict[str, Any]:
    if not isinstance(config, dict):
        raise ValueError("config must be a dictionary")
    unknown = set(config) - CONFIG_KEYS
    if unknown:
        raise ValueError(f"Unknown config keys: {', '.join(sorted(map(str, unknown)))}")
    for key in ("reference_gene", "calibrator_group", "calibrator_time_h"):
        if key not in config:
            raise ValueError(f"Missing required config key: {key}")
    result = {
        "reference_gene": _identifier(config["reference_gene"], "reference_gene"),
        "calibrator_group": _identifier(config["calibrator_group"], "calibrator_group"),
        "calibrator_time_h": _finite_number(config["calibrator_time_h"], "calibrator_time_h", minimum=0),
        "period_h": _finite_number(config.get("period_h", 24), "period_h", positive=True),
        "design": config.get("design", "unknown"),
        "technical_sd_limit": _finite_number(config.get("technical_sd_limit", 0.5), "technical_sd_limit", minimum=0),
    }
    if result["design"] not in ("unknown", "independent", "repeated"):
        raise ValueError("design must be one of: unknown, independent, repeated")
    return result


def _qc(qc: list[dict[str, Any]], level: str, code: str, message: str, **context: Any) -> None:
    qc.append({"level": level, "code": code, "message": message, **context})


def _mean(values: list[float]) -> float:
    # Scale first to avoid overflowing the sum or underflowing every summand.
    scale = max(abs(value) for value in values)
    if scale == 0:
        return 0.0
    return scale * (math.fsum(value / scale for value in values) / len(values))


def _sample_sd(values: list[float], mean: float) -> float | None:
    if len(values) < 2:
        return None
    scale = max(abs(value) for value in values)
    if scale == 0:
        return 0.0
    variance_scaled = math.fsum((value / scale - mean / scale) ** 2 for value in values) / (len(values) - 1)
    result = scale * math.sqrt(variance_scaled)
    return result if math.isfinite(result) else None


def _read_csv(csv_path: str | Path) -> tuple[bytes, list[dict[str, str]], list[dict[str, Any]], dict[str, dict[str, Any]], list[str]]:
    payload = Path(csv_path).read_bytes()
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError("CSV must use UTF-8 (an optional UTF-8 BOM is supported)") from exc
    try:
        rows = list(csv.reader(io.StringIO(text, newline=""), strict=True))
    except csv.Error as exc:
        raise ValueError(f"Malformed CSV: {exc}") from exc
    if not rows:
        raise ValueError("CSV is empty")
    header = rows[0]
    if any(not column for column in header):
        raise ValueError("CSV contains an empty column name")
    if len(header) != len(set(header)):
        raise ValueError("CSV column names must be unique")
    missing = set(REQUIRED_COLUMNS) - set(header)
    if missing:
        raise ValueError(f"Missing required CSV columns: {', '.join(sorted(missing))}")
    raw_rows: list[dict[str, str]] = []
    parsed: list[dict[str, Any]] = []
    samples: dict[str, dict[str, Any]] = {}
    genes: list[str] = []
    well_keys: set[tuple[str, str, str]] = set()
    biological_keys: dict[tuple[str, str, float], str] = {}
    subject_groups: dict[str, str] = {}
    for row_number, cells in enumerate(rows[1:], start=2):
        if not cells:  # Conventional empty CSV lines are not observations.
            continue
        if len(cells) != len(header):
            raise ValueError(f"CSV row {row_number} has {len(cells)} fields; expected {len(header)}")
        raw = dict(zip(header, cells))
        row: dict[str, Any] = dict(raw)
        for name in ("sample_id", "subject_id", "group", "gene", "tech_rep", "run_id"):
            row[name] = _identifier(raw[name], f"row {row_number} {name}")
        row["time_h"] = _finite_number(raw["time_h"], f"row {row_number} time_h", minimum=0)
        ct_text = raw["ct"].strip()
        row["ct"] = None if ct_text.casefold() in MISSING_CT_TOKENS else _finite_number(ct_text, f"row {row_number} ct")
        row["source_row"] = row_number
        well_key = (row["sample_id"], row["gene"], row["tech_rep"])
        if well_key in well_keys:
            raise ValueError(f"Duplicate sample_id/gene/tech_rep at CSV row {row_number}: {well_key!r}")
        well_keys.add(well_key)
        metadata = {name: row[name] for name in ("sample_id", "subject_id", "group", "time_h", "run_id")}
        if row["sample_id"] in samples and samples[row["sample_id"]] != metadata:
            raise ValueError(f"Inconsistent metadata (including run_id) for sample_id {row['sample_id']!r}")
        samples[row["sample_id"]] = metadata
        bio_key = (row["subject_id"], row["group"], row["time_h"])
        prior_sample = biological_keys.get(bio_key)
        if prior_sample is not None and prior_sample != row["sample_id"]:
            raise ValueError(f"Multiple sample_id values for the same subject/group/time: {bio_key!r}; combine technical material before analysis")
        biological_keys[bio_key] = row["sample_id"]
        prior_group = subject_groups.get(row["subject_id"])
        if prior_group is not None and prior_group != row["group"]:
            raise ValueError(f"subject_id {row['subject_id']!r} occurs in more than one group")
        subject_groups[row["subject_id"]] = row["group"]
        if row["gene"] not in genes:
            genes.append(row["gene"])
        raw_rows.append(raw)
        parsed.append(row)
    if not parsed:
        raise ValueError("CSV has no data rows")
    return payload, raw_rows, parsed, samples, genes


def _phase_count(times: list[float], period: float) -> int:
    fractions = sorted((time % period) / period for time in times)
    distinct: list[float] = []
    for fraction in fractions:
        if not distinct or fraction - distinct[-1] > 1e-8:
            distinct.append(fraction)
    if len(distinct) > 1 and distinct[0] + 1 - distinct[-1] <= 1e-8:
        distinct.pop()
    return len(distinct)


def _fit_group(rows: list[dict[str, Any]], gene: str, group: str,
               config: dict[str, Any], global_reasons: list[str],
               qc: list[dict[str, Any]]) -> dict[str, Any]:
    valid = [row for row in rows if row["log2rq"] is not None]
    period = config["period_h"]
    times = [row["time_h"] for row in valid]
    distinct_times = sorted(set(times))
    input_times = sorted({row["time_h"] for row in rows})
    biological_counts = [
        {"time_h": time, "n_subjects": len({row["subject_id"] for row in valid if row["time_h"] == time})}
        for time in input_times
    ]
    fit: dict[str, Any] = {
        "gene": gene, "group": group, "status": "not_testable", "reason": None,
        "n": len(valid), "n_input": len(rows), "n_valid": len(valid),
        "n_subjects": len({row["subject_id"] for row in valid}),
        "n_times": len(distinct_times), "n_times_input": len(input_times),
        "n_phases": _phase_count(times, period),
        "mesor": None, "cos_coef": None, "sin_coef": None, "amplitude": None,
        "peak_h": None, "p_value": None, "q_value": None, "period_h": period,
        "r_squared": None, "residual_df": len(valid) - 3,
        "f_statistic": None, "condition_number": None,
        "biological_replicates_by_time": biological_counts,
    }
    excluded = len(rows) - len(valid)
    if excluded:
        _qc(qc, "warning", "missing_observations_excluded",
            "Incomplete normalized observations were excluded from the fit; missingness may be nonrandom and is not evidence of no change.",
            gene=gene, group=group, n_input=len(rows), n_valid=len(valid), n_excluded=excluded)
    if not valid:
        fit["reason"] = "no_complete_observations"
        return fit
    phase = np.asarray([(time % period) / period for time in times], dtype=float) * (2 * math.pi)
    design_matrix = np.column_stack((np.ones(len(valid)), np.cos(phase), np.sin(phase)))
    if np.linalg.matrix_rank(design_matrix) < 3:
        fit["reason"] = "rank_deficient_design"
        return fit
    condition = float(np.linalg.cond(design_matrix))
    fit["condition_number"] = condition if math.isfinite(condition) else None
    if not math.isfinite(condition) or condition > 1e8:
        fit["reason"] = "ill_conditioned_design"
        return fit
    y = np.asarray([row["log2rq"] for row in valid], dtype=float)
    scale = max(1.0, float(np.max(np.abs(y))))
    y_scaled = y / scale
    beta_scaled = np.linalg.lstsq(design_matrix, y_scaled, rcond=None)[0]
    with np.errstate(over="ignore", invalid="ignore"):
        beta = beta_scaled * scale
    if not np.all(np.isfinite(beta)):
        fit["reason"] = "fit_numeric_range_error"
        return fit
    mesor, cosine, sine = map(float, beta)
    amplitude = math.hypot(cosine, sine)
    if not math.isfinite(amplitude):
        fit["reason"] = "fit_numeric_range_error"
        return fit
    constant_response = float(np.max(y_scaled) - np.min(y_scaled)) <= 128 * np.finfo(float).eps
    negligible_amplitude = amplitude / scale <= 128 * np.finfo(float).eps
    if constant_response:
        mesor, cosine, sine, amplitude = _mean(y.tolist()), 0.0, 0.0, 0.0
    elif negligible_amplitude:
        amplitude = 0.0
    peak = None if amplitude == 0 else ((math.atan2(sine, cosine) % (2 * math.pi)) / (2 * math.pi)) * period
    if peak is not None and math.isclose(peak, period, rel_tol=0, abs_tol=period * 1e-10):
        peak = 0.0
    residual = y_scaled - design_matrix @ beta_scaled
    sse = float(residual @ residual)
    centered = y_scaled - float(np.mean(y_scaled))
    sst = float(centered @ centered)
    fit.update(mesor=mesor, cos_coef=cosine, sin_coef=sine, amplitude=amplitude,
               peak_h=peak, r_squared=None if constant_response else max(0.0, min(1.0, 1.0 - sse / sst)))
    if fit["residual_df"] <= 0:
        fit["reason"] = "insufficient_residual_degrees_of_freedom"
        return fit
    reasons = list(global_reasons)
    if config["design"] != "independent":
        reasons.insert(0, f"design_{config['design']}_descriptive_only")
    if fit["n_phases"] < 4:
        reasons.append("fewer_than_four_distinct_phases")
    if max(times) - min(times) < 0.75 * period:
        reasons.append("time_span_below_three_quarters_period")
    if any(item["n_subjects"] < 2 for item in biological_counts):
        reasons.append("fewer_than_two_subjects_at_a_time")
    if any(not row["inference_eligible"] for row in valid):
        reasons.append("normalization_quality_gate_failed")
    perfect_fit = not constant_response and sse <= (128 * np.finfo(float).eps) ** 2 * max(1.0, sst)
    if perfect_fit:
        reasons.append("perfect_fit_no_residual_variance")
    if reasons:
        fit["status"] = "descriptive"
        fit["reason"] = "; ".join(dict.fromkeys(reasons))
        _qc(qc, "info", "inferential_test_withheld",
            "Cosinor estimates are descriptive; no P value or Q value is reported for this fit.",
            gene=gene, group=group, reason=fit["reason"], n_input=len(rows), n_valid=len(valid), n_subjects=fit["n_subjects"])
        return fit
    if constant_response:
        f_statistic, p_value = 0.0, 1.0
    else:
        model_ss = max(0.0, sst - sse)
        f_statistic = (model_ss / 2.0) / (sse / fit["residual_df"])
        p_value = float(f_distribution.sf(f_statistic, 2, fit["residual_df"]))
    if not math.isfinite(f_statistic) or not math.isfinite(p_value):
        fit["status"], fit["reason"] = "descriptive", "test_numeric_range_error"
        return fit
    fit.update(status="tested", reason=None, f_statistic=float(f_statistic), p_value=p_value)
    return fit


def _adjust_bh(fits: list[dict[str, Any]]) -> int:
    tested = [(index, fit["p_value"]) for index, fit in enumerate(fits) if fit["p_value"] is not None]
    tested.sort(key=lambda entry: entry[1])
    running = 1.0
    for rank in range(len(tested), 0, -1):
        index, p_value = tested[rank - 1]
        running = min(running, p_value * len(tested) / rank)
        fits[index]["q_value"] = float(min(1.0, running))
    return len(tested)


def analyze(csv_path: str | Path, config: dict[str, Any]) -> dict[str, Any]:
    """Analyze a UTF-8 long-form Ct CSV without changing source identifiers.

    Structural errors raise ``ValueError``. Expected missing measurements stay
    explicit in ``normalized`` with ``None`` values and reasons. ``raw_rows``
    contains the original CSV strings, including unrecognized extra columns.
    ``n`` in a fit counts sample-level normalized observations, never wells.
    ``inference_eligible`` on normalized rows is a technical/normalization gate;
    the design, sampling and single-run gates are additionally enforced on fits.
    """
    cfg = _validate_config(config)
    payload, raw_rows, parsed, samples, genes = _read_csv(csv_path)
    qc: list[dict[str, Any]] = []
    reference = cfg["reference_gene"]
    targets = [gene for gene in genes if gene != reference]
    groups = list(dict.fromkeys(sample["group"] for sample in samples.values()))
    runs = list(dict.fromkeys(row["run_id"] for row in parsed))
    missing_ct_count = sum(row["ct"] is None for row in parsed)
    if missing_ct_count:
        _qc(qc, "warning", "missing_ct",
            "Missing/Undetermined Ct values remain missing and are never replaced with zero; missingness may be nonrandom.",
            n_missing=missing_ct_count, n_raw_rows=len(parsed))
    _qc(qc, "info", "normalization_assumptions",
        "Normalization assumes one stable reference gene and equal amplification efficiencies; Ct values alone do not validate these assumptions.")
    if reference not in genes:
        _qc(qc, "warning", "reference_gene_absent", "The configured reference gene is absent; normalized values cannot be calculated.", reference_gene=reference)
    if not targets:
        _qc(qc, "warning", "no_target_genes", "No target genes are present after excluding the reference gene.")
    well_groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in parsed:
        well_groups[(row["sample_id"], row["gene"])].append(row)
    aggregates: dict[tuple[str, str], dict[str, Any]] = {}
    for (sample_id, gene), wells in well_groups.items():
        values = [well["ct"] for well in wells if well["ct"] is not None]
        mean = _mean(values) if values else None
        sd = _sample_sd(values, mean) if mean is not None else None
        sd_range_error = len(values) >= 2 and sd is None
        high_sd = sd is not None and sd > cfg["technical_sd_limit"]
        aggregates[(sample_id, gene)] = {
            "mean_ct": mean, "sd_ct": sd, "n_tech": len(values), "n_tech_input": len(wells),
            "quality_failed": bool(high_sd or sd_range_error),
        }
        if high_sd:
            _qc(qc, "warning", "high_technical_sd", "Technical-replicate Ct SD exceeds the configured limit; values are retained for description but related inferential tests are withheld.",
                sample_id=sample_id, gene=gene, sd_ct=sd, technical_sd_limit=cfg["technical_sd_limit"], n_tech=len(values))
        elif sd_range_error:
            _qc(qc, "warning", "technical_sd_numeric_range_error", "Technical-replicate SD is outside the finite numeric range; related inferential tests are withheld.", sample_id=sample_id, gene=gene)
        if len(values) == 1:
            _qc(qc, "warning", "single_valid_technical_replicate", "Only one valid technical replicate remains; technical SD cannot be estimated and the well is not an additional biological replicate.", sample_id=sample_id, gene=gene)
    normalized: list[dict[str, Any]] = []
    empty = {"mean_ct": None, "sd_ct": None, "n_tech": 0, "n_tech_input": 0, "quality_failed": False}
    for gene in targets:
        gene_rows: list[dict[str, Any]] = []
        for sample_id, metadata in samples.items():
            target = aggregates.get((sample_id, gene), empty)
            ref = aggregates.get((sample_id, reference), empty)
            reasons: list[str] = []
            if target["mean_ct"] is None:
                reasons.append("missing_target_ct")
            if ref["mean_ct"] is None:
                reasons.append("missing_reference_ct")
            delta = None if reasons else target["mean_ct"] - ref["mean_ct"]
            if delta is not None and not math.isfinite(delta):
                delta = None
                reasons.append("delta_ct_numeric_range_error")
            quality_reasons = []
            if target["quality_failed"]:
                quality_reasons.append("target_technical_qc_failed")
            if ref["quality_failed"]:
                quality_reasons.append("reference_technical_qc_failed")
            row = {
                **metadata, "gene": gene, "target_mean_ct": target["mean_ct"],
                "reference_mean_ct": ref["mean_ct"], "target_sd_ct": target["sd_ct"],
                "reference_sd_ct": ref["sd_ct"], "target_n_tech": target["n_tech"],
                "reference_n_tech": ref["n_tech"], "target_n_tech_input": target["n_tech_input"],
                "reference_n_tech_input": ref["n_tech_input"], "delta_ct": delta,
                "calibrator_delta_ct": None, "calibrator_n": 0, "log2rq": None, "rq": None,
                "status": "missing", "reason": "; ".join(reasons) if reasons else None,
                "inference_eligible": delta is not None and not target["quality_failed"] and not ref["quality_failed"],
                "inference_reason": "; ".join(quality_reasons) if quality_reasons else None,
            }
            gene_rows.append(row)
        calibrator_candidates = [row for row in gene_rows if row["group"] == cfg["calibrator_group"] and row["time_h"] == cfg["calibrator_time_h"]]
        calibrators = [row for row in calibrator_candidates if row["delta_ct"] is not None]
        calibrator_mean = _mean([row["delta_ct"] for row in calibrators]) if calibrators else None
        calibrator_qc_pass = all(row["inference_eligible"] for row in calibrators)
        if not calibrators:
            _qc(qc, "warning", "missing_calibrator", "No complete target/reference pair is available at the configured calibrator group and time.", gene=gene, calibrator_group=cfg["calibrator_group"], calibrator_time_h=cfg["calibrator_time_h"], n_candidates=len(calibrator_candidates))
        elif len(calibrators) < len(calibrator_candidates):
            _qc(qc, "warning", "incomplete_calibrators_excluded", "Incomplete calibrator observations were excluded from the mean; calibrator missingness remains visible.", gene=gene, n_input=len(calibrator_candidates), n_valid=len(calibrators))
        for row in gene_rows:
            row["calibrator_delta_ct"] = calibrator_mean
            row["calibrator_n"] = len(calibrators)
            reasons = [row["reason"]] if row["reason"] else []
            if calibrator_mean is None:
                reasons.append("missing_calibrator")
                row["inference_eligible"] = False
            elif row["delta_ct"] is not None:
                log2rq = -(row["delta_ct"] - calibrator_mean)
                if not math.isfinite(log2rq):
                    reasons.append("log2rq_numeric_range_error")
                    row["inference_eligible"] = False
                else:
                    row["log2rq"] = log2rq
                    try:
                        rq = 2.0 ** log2rq
                    except OverflowError:
                        rq = math.inf
                    if not math.isfinite(rq) or rq <= 0:
                        row["status"] = "partial"
                        reasons.append("rq_out_of_range")
                        row["inference_eligible"] = False
                        _qc(qc, "warning", "rq_out_of_range", "RQ is outside the finite positive numeric range; log2RQ is retained for description and inferential testing is withheld.", sample_id=row["sample_id"], gene=gene)
                    else:
                        row["rq"] = rq
                        row["status"] = "ok"
                    row["inference_eligible"] = row["inference_eligible"] and calibrator_qc_pass
            row["reason"] = "; ".join(reasons) if reasons else None
            inference_reasons = [reason for reason in (row["inference_reason"], row["reason"]) if reason]
            if calibrators and not calibrator_qc_pass:
                inference_reasons.append("calibrator_technical_qc_failed")
            row["inference_reason"] = "; ".join(inference_reasons) if inference_reasons else None
            normalized.append(row)
    subject_times: dict[str, set[float]] = defaultdict(set)
    for sample in samples.values():
        subject_times[sample["subject_id"]].add(sample["time_h"])
    global_reasons: list[str] = []
    if len(runs) != 1:
        global_reasons.append("multiple_runs_without_validated_batch_model")
        _qc(qc, "warning", "multiple_runs", "More than one run_id is present; no inter-run calibration or batch model is implemented, so all fits are descriptive.", run_ids=runs)
    repeated_subjects = [subject for subject, times in subject_times.items() if len(times) > 1]
    if repeated_subjects:
        global_reasons.append("subject_repeated_across_time")
        _qc(qc, "warning", "repeated_subjects", "Some subjects occur at multiple times; ordinary independent-error F tests are withheld for this dataset.", n_repeated_subjects=len(repeated_subjects))
    fits = [_fit_group([row for row in normalized if row["gene"] == gene and row["group"] == group], gene, group, cfg, global_reasons, qc)
            for gene in targets for group in groups]
    family_size = _adjust_bh(fits)
    _qc(qc, "info", "sample_level_analysis", "Technical wells were averaged within sample and gene before normalization; fit n counts sample-level observations, and n_subjects counts distinct biological subjects.", n_samples=len(samples), n_subjects=len(subject_times))
    return {
        "raw_rows": raw_rows, "normalized": normalized, "fits": fits, "qc": qc,
        "config": cfg, "source_sha256": hashlib.sha256(payload).hexdigest(),
        "audit": {
            "required_columns": list(REQUIRED_COLUMNS), "n_raw_rows": len(raw_rows),
            "n_samples": len(samples), "n_subjects": len(subject_times), "n_missing_ct": missing_ct_count,
            "missing_ct_tokens_case_insensitive": sorted(MISSING_CT_TOKENS),
            "identifier_policy": "Identifiers preserve original case and whitespace; only empty or whitespace-only identifiers are rejected.",
            "nonfinite_policy": "NaN and infinity are invalid input, not missing Ct tokens.",
            "technical_aggregation": "Arithmetic mean Ct of available technical wells; sample SD uses n-1.",
            "normalization": "delta_ct = target mean Ct - reference mean Ct; log2rq = -(delta_ct - mean complete calibrator delta_ct); rq = 2**log2rq.",
            "calibrator_policy": "Mean of complete sample-level delta Ct at exactly the configured group and time; no imputation.",
            "fit_unit": "One sample-level normalized observation; technical wells never increase biological N.",
            "normalized_inference_eligible_definition": "Technical/normalization quality gate only; design and sampling gates are additional fit-level requirements.",
            "model": "log2rq = mesor + cos_coef*cos(2*pi*time_h/period_h) + sin_coef*sin(2*pi*time_h/period_h)",
            "peak_origin": "Hours after input time_h=0, modulo period_h; no peak is reported at numerically zero amplitude.",
            "inference_requirements": ["design=independent", "one run_id in the entire input", "no subject sampled at multiple times", "at least two valid distinct subjects at every input time for that group, including times with wholly missing measurements", "at least four distinct phases modulo period", "observed time span >= 0.75*period", "full-rank design with condition number <= 1e8", "positive residual degrees of freedom", "no related technical/normalization quality failure", "finite residual variance for nonconstant response"],
            "bh_family": "All finite zero-amplitude F-test P values across every target gene/group fit in this analyze invocation.",
            "bh_family_size": family_size, "run_ids": runs,
            "n_repeated_subjects": len(repeated_subjects), "phase_tolerance_fraction_of_period": 1e-8,
            "group_difference_tested": False,
        },
    }
