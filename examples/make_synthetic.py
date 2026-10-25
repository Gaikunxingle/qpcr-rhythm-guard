"""Deterministic artificial values; no real sample or patient records."""
import csv
import math
import random
from pathlib import Path


def main():
    rng = random.Random(20261025)
    dest = Path(__file__).with_name("synthetic.csv")
    fields = ["sample_id", "subject_id", "group", "time_h", "gene", "tech_rep", "ct", "run_id"]
    rows = []
    for group in ("Control", "Perturbed"):
        for time in range(0, 24, 4):
            for bio in range(1, 4):
                sample = f"SYN_{group}_{time:02d}_{bio}"
                ref = 19 + rng.gauss(0, 0.2)
                expression = (1.2 if group == "Control" else 0.5) * math.cos(2 * math.pi * (time - 6) / 24)
                # StableB has balanced offsets at each time, serving as a synthetic flat control.
                levels = {"GAPDH": ref, "ClockA": ref + 5 - expression + rng.gauss(0, 0.14), "StableB": ref + 6 + .2 * (bio - 2) + rng.gauss(0, 0.02)}
                for gene, level in levels.items():
                    for tech in (1, 2):
                        ct = f"{level + rng.gauss(0, 0.05):.6f}"
                        if group == "Perturbed" and time == 12 and bio == 3 and gene == "StableB":
                            ct = "Undetermined"
                        rows.append(dict(zip(fields, [sample, sample, group, time, gene, tech, ct, "SYN_RUN_1"])))
    with dest.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Created {len(rows)} artificial well records: {dest}")


if __name__ == "__main__":
    main()
