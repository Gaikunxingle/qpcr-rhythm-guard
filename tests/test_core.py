"""Independent numerical oracles and scientific guardrails for the public API.

The principal oracle has six equally spaced phases, two independent biological
samples per phase, and paired residuals of +/-0.2. Orthogonality gives an exact
signal sum of squares and residual sum of squares, without reproducing the
implementation's least-squares procedure.
"""

import csv
import math
from pathlib import Path
import tempfile
import unittest

from scipy.stats import f as f_distribution

from qpcr_rhythm_guard.core import analyze


FIELDS = [
    "sample_id", "subject_id", "group", "time_h", "gene", "tech_rep",
    "ct", "run_id",
]


def synthetic_rows(amplitudes=None, times=None, biological_reps=2,
                   technical_reps=2, repeated=False):
    """Build Ct-level observations with varying reference Ct and known signal."""
    amplitudes = amplitudes or {"Clock": 1.5}
    times = [0, 4, 8, 12, 16, 20] if times is None else times
    rows = []
    for phase_index, time_h in enumerate(times):
        for biological_rep in range(biological_reps):
            sample_id = f"sample_{phase_index}_{biological_rep}"
            subject_id = (f"subject_{biological_rep}" if repeated else sample_id)
            # Varying reference expression tests sample-level delta-Ct pairing.
            reference_ct = 19.0 + 0.1 * phase_index + 0.03 * biological_rep
            residual = (0.2 * (2 * biological_rep - biological_reps + 1)
                        if biological_reps > 1 else 0.0)
            for technical_rep in range(technical_reps):
                technical_offset = 0.02 * (2 * technical_rep - technical_reps + 1)
                base = dict(sample_id=sample_id, subject_id=subject_id,
                            group="Control", time_h=time_h,
                            tech_rep=technical_rep + 1, run_id="run_1")
                rows.append(dict(base, gene="Reference",
                                 ct=reference_ct + technical_offset))
                for gene, amplitude in amplitudes.items():
                    signal = amplitude * math.cos(2 * math.pi * (time_h - 4) / 24)
                    target_ct = reference_ct + 5 - signal - residual
                    rows.append(dict(base, gene=gene,
                                     ct=target_ct + technical_offset))
    return rows


class CoreTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.input_path = Path(self.temporary_directory.name) / "input.csv"
        self.config = {
            "reference_gene": "Reference",
            "calibrator_group": "Control",
            "calibrator_time_h": 0,
            "period_h": 24,
            "design": "independent",
            "technical_sd_limit": 0.5,
        }

    def run_analysis(self, rows=None, **overrides):
        rows = synthetic_rows() if rows is None else rows
        with self.input_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerows(rows)
        return analyze(self.input_path, dict(self.config, **overrides))

    def test_sample_normalization_and_known_amplitude_peak(self):
        result = self.run_analysis()
        normalized = result["normalized"]
        self.assertEqual(len(normalized), 12)
        for row in normalized:
            biological_rep = int(row["sample_id"].split("_")[-1])
            residual = -0.2 if biological_rep == 0 else 0.2
            signal = 1.5 * math.cos(2 * math.pi * (row["time_h"] - 4) / 24)
            expected_log2rq = signal + residual - 0.75
            self.assertEqual(row["status"], "ok")
            self.assertAlmostEqual(row["delta_ct"], 5 - signal - residual, places=11)
            self.assertAlmostEqual(row["calibrator_delta_ct"], 4.25, places=11)
            self.assertAlmostEqual(row["log2rq"], expected_log2rq, places=11)
            self.assertAlmostEqual(row["rq"], 2 ** expected_log2rq, places=11)
        fit = result["fits"][0]
        self.assertEqual(fit["status"], "tested")
        self.assertEqual(fit["n"], 12)
        self.assertEqual(fit["n_subjects"], 12)
        self.assertAlmostEqual(fit["mesor"], -0.75, places=11)
        self.assertAlmostEqual(fit["amplitude"], 1.5, places=11)
        self.assertAlmostEqual(fit["peak_h"], 4.0, places=11)

    def test_joint_model_test_matches_orthogonal_signal_oracle(self):
        fit = self.run_analysis()["fits"][0]
        # N=12: signal SS=N*A^2/2=13.5; error SS=N*0.2^2=0.48.
        # The sinusoidal model adds two coefficients to an estimated intercept.
        expected_f = (13.5 / 2) / (0.48 / 9)
        expected_p = float(f_distribution.sf(expected_f, 2, 9))
        self.assertAlmostEqual(fit["p_value"], expected_p, delta=expected_p * 1e-10)
        self.assertAlmostEqual(fit["q_value"], expected_p, delta=expected_p * 1e-10)
        self.assertAlmostEqual(fit["r_squared"], 13.5 / 13.98, places=11)

    def test_additional_technical_wells_do_not_increase_biological_n(self):
        two = self.run_analysis(synthetic_rows(technical_reps=2))["fits"][0]
        four = self.run_analysis(synthetic_rows(technical_reps=4))["fits"][0]
        self.assertEqual(two["n"], 12)
        self.assertEqual(four["n"], 12)
        self.assertEqual(four["n_subjects"], 12)
        self.assertAlmostEqual(two["amplitude"], four["amplitude"], places=11)
        self.assertAlmostEqual(two["p_value"], four["p_value"], places=14)

    def test_technical_wells_cannot_substitute_for_biological_replication(self):
        fit = self.run_analysis(synthetic_rows(biological_reps=1,
                                               technical_reps=4))["fits"][0]
        self.assertEqual(fit["n"], 6)
        self.assertEqual(fit["status"], "descriptive")
        self.assertIsNone(fit["p_value"])
        self.assertIsNone(fit["q_value"])

    def test_unknown_design_remains_descriptive(self):
        fit = self.run_analysis(design="unknown")["fits"][0]
        self.assertEqual(fit["status"], "descriptive")
        self.assertAlmostEqual(fit["amplitude"], 1.5, places=11)
        self.assertIsNone(fit["p_value"])
        self.assertIsNone(fit["q_value"])

    def test_constant_response_has_no_artificial_peak(self):
        rows = synthetic_rows()
        references = {(row["sample_id"], row["tech_rep"]): row["ct"]
                      for row in rows if row["gene"] == "Reference"}
        for row in rows:
            if row["gene"] == "Clock":
                row["ct"] = references[(row["sample_id"], row["tech_rep"])] + 5
        fit = self.run_analysis(rows)["fits"][0]
        self.assertEqual(fit["status"], "tested")
        self.assertEqual(fit["amplitude"], 0)
        self.assertIsNone(fit["peak_h"])
        self.assertEqual(fit["p_value"], 1)
        self.assertEqual(fit["q_value"], 1)

    def test_perfect_nonconstant_signal_withholds_zero_residual_variance_test(self):
        rows = synthetic_rows()
        for row in rows:
            if row["gene"] == "Clock":
                biological_rep = int(row["sample_id"].split("_")[-1])
                row["ct"] += -0.2 if biological_rep == 0 else 0.2
        fit = self.run_analysis(rows)["fits"][0]
        self.assertEqual(fit["status"], "descriptive")
        self.assertAlmostEqual(fit["amplitude"], 1.5, places=11)
        self.assertAlmostEqual(fit["peak_h"], 4.0, places=11)
        self.assertIn("perfect_fit_no_residual_variance", fit["reason"])
        self.assertIsNone(fit["p_value"])
        self.assertIsNone(fit["q_value"])

    def test_repeated_subjects_do_not_receive_independent_model_p_values(self):
        for design in ["repeated", "independent"]:
            with self.subTest(design=design):
                fit = self.run_analysis(synthetic_rows(repeated=True),
                                        design=design)["fits"][0]
                self.assertEqual(fit["status"], "descriptive")
                self.assertEqual(fit["n_subjects"], 2)
                self.assertIsNone(fit["p_value"])
                self.assertIsNone(fit["q_value"])

    def test_multiple_sample_ids_from_one_subject_are_not_independent(self):
        rows = synthetic_rows()
        # Two aliquots from each subject are distinct samples, not new subjects.
        additional_samples = [dict(row, sample_id=row["sample_id"] + "_aliquot")
                              for row in rows]
        with self.assertRaises(ValueError):
            self.run_analysis(rows + additional_samples)

    def test_missing_target_and_reference_are_preserved_not_zero_filled(self):
        rows = synthetic_rows()
        for row in rows:
            if row["sample_id"] == "sample_2_0" and row["gene"] == "Clock":
                row["ct"] = "Undetermined"
            if row["sample_id"] == "sample_2_1" and row["gene"] == "Reference":
                row["ct"] = ""
        result = self.run_analysis(rows)
        missing = [row for row in result["normalized"]
                   if row["sample_id"] in {"sample_2_0", "sample_2_1"}]
        self.assertEqual(len(missing), 2)
        for row in missing:
            self.assertEqual(row["status"], "missing")
            self.assertIsNone(row["log2rq"])
            self.assertIsNone(row["rq"])
        self.assertEqual(result["fits"][0]["n"], 10)
        self.assertTrue(result["qc"])

    def test_one_missing_technical_well_uses_the_remaining_well(self):
        rows = synthetic_rows()
        for row in rows:
            if (row["sample_id"] == "sample_2_0" and row["gene"] == "Clock"
                    and row["tech_rep"] == 1):
                row["ct"] = "Undetermined"
        normalized = self.run_analysis(rows)["normalized"]
        observed = next(row for row in normalized if row["sample_id"] == "sample_2_0")
        self.assertEqual(observed["status"], "ok")
        self.assertEqual(observed["target_n_tech"], 1)
        self.assertEqual(observed["reference_n_tech"], 2)
        # Remaining target well is +0.02 Ct: expression shifts down by 0.02.
        self.assertAlmostEqual(observed["log2rq"], -0.22, places=11)

    def test_missing_calibrator_does_not_fall_back_to_other_samples(self):
        rows = synthetic_rows()
        for row in rows:
            if row["time_h"] == 0 and row["gene"] == "Clock":
                row["ct"] = "Undetermined"
        result = self.run_analysis(rows)
        self.assertTrue(result["normalized"])
        self.assertTrue(all(row["log2rq"] is None for row in result["normalized"]))
        self.assertEqual(result["fits"][0]["status"], "not_testable")
        self.assertIsNone(result["fits"][0]["p_value"])

    def test_completely_missing_time_point_remains_in_replication_gate(self):
        rows = synthetic_rows()
        for row in rows:
            if row["time_h"] == 8 and row["gene"] == "Clock":
                row["ct"] = "Undetermined"
        fit = self.run_analysis(rows)["fits"][0]
        self.assertEqual(fit["n"], 10)
        self.assertEqual(fit["status"], "descriptive")
        counts = {row["time_h"]: row["n_subjects"]
                  for row in fit["biological_replicates_by_time"]}
        self.assertEqual(counts[8], 0)
        self.assertIsNone(fit["p_value"])
        self.assertIsNone(fit["q_value"])

    def test_high_technical_dispersion_blocks_inference(self):
        rows = synthetic_rows()
        for row in rows:
            if row["sample_id"] == "sample_2_0" and row["gene"] == "Clock":
                row["ct"] += -1.0 if row["tech_rep"] == 1 else 1.0
        result = self.run_analysis(rows)
        self.assertEqual(result["fits"][0]["status"], "descriptive")
        self.assertIsNone(result["fits"][0]["p_value"])
        self.assertTrue(result["qc"])

    def test_unmodelled_multiple_runs_block_inference(self):
        rows = synthetic_rows()
        for row in rows:
            if row["time_h"] >= 12:
                row["run_id"] = "run_2"
        fit = self.run_analysis(rows)["fits"][0]
        self.assertEqual(fit["status"], "descriptive")
        self.assertIsNone(fit["p_value"])

    def test_short_time_span_is_descriptive(self):
        fit = self.run_analysis(synthetic_rows(times=[0, 4, 8, 12]))["fits"][0]
        self.assertEqual(fit["status"], "descriptive")
        self.assertIsNone(fit["p_value"])

    def test_singular_phase_design_cannot_report_a_fitted_rhythm(self):
        # Four different elapsed times, but only two phases modulo 24 hours.
        fit = self.run_analysis(synthetic_rows(times=[0, 12, 24, 36]))["fits"][0]
        self.assertEqual(fit["status"], "not_testable")
        self.assertIsNone(fit["amplitude"])
        self.assertIsNone(fit["peak_h"])
        self.assertIsNone(fit["p_value"])

    def test_bh_family_matches_known_four_hypothesis_oracle(self):
        desired_p = {"ClockA": .01, "ClockB": .04, "ClockC": .03, "ClockD": .20}
        # Invert the known F oracle to construct four signals with chosen P.
        amplitudes = {gene: math.sqrt(float(f_distribution.isf(p, 2, 9)) * .16 / 9)
                      for gene, p in desired_p.items()}
        fits = {fit["gene"]: fit
                for fit in self.run_analysis(synthetic_rows(amplitudes))["fits"]}
        expected_q = {"ClockA": .04, "ClockB": .16 / 3,
                      "ClockC": .16 / 3, "ClockD": .20}
        self.assertEqual(set(fits), set(desired_p))
        for gene in desired_p:
            self.assertAlmostEqual(fits[gene]["p_value"], desired_p[gene], places=10)
            self.assertAlmostEqual(fits[gene]["q_value"], expected_q[gene], places=10)

    def test_duplicate_well_identity_is_rejected(self):
        rows = synthetic_rows()
        rows.append(dict(rows[0]))
        with self.assertRaises(ValueError):
            self.run_analysis(rows)

    def test_nonfinite_numeric_input_is_rejected(self):
        for field in ["ct", "time_h"]:
            for value in ["inf", "-inf", "nan"]:
                with self.subTest(field=field, value=value):
                    rows = synthetic_rows()
                    rows[0][field] = value
                    with self.assertRaises(ValueError):
                        self.run_analysis(rows)

    def test_unknown_configuration_key_is_not_silently_ignored(self):
        with self.assertRaises(ValueError):
            self.run_analysis(technical_sd_limt=0.25)

    def test_conflicting_metadata_for_one_sample_is_rejected(self):
        rows = synthetic_rows()
        rows[0]["subject_id"] = "different_subject"
        with self.assertRaises(ValueError):
            self.run_analysis(rows)


if __name__ == "__main__":
    unittest.main()
