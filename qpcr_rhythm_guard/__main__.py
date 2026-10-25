"""Local CLI; never uploads data or calls an AI API."""
import argparse
import json
import sys
from pathlib import Path

from .core import analyze
from .reporting import export_results


def main():
    parser = argparse.ArgumentParser(description="qPCR QC and fixed-period cosinor; local processing")
    parser.add_argument("input", type=Path, help="UTF-8 long-table CSV")
    parser.add_argument("--config", type=Path, required=True, help="JSON analysis configuration")
    parser.add_argument("--out", type=Path, required=True, help="New or empty output directory")
    args = parser.parse_args()
    try:
        if args.out.exists() and (not args.out.is_dir() or any(args.out.iterdir())):
            raise ValueError("Output must be a new or empty directory; existing results are preserved.")
        with args.config.open(encoding="utf-8-sig") as stream:
            config = json.load(stream)
        if not isinstance(config, dict):
            raise ValueError("Configuration must be a JSON object.")
        result = analyze(args.input, config)
        export_results(result, args.out, args.input, args.config)
    except (ValueError, OSError, KeyError) as exc:
        print(f"Cannot analyze: {exc}", file=sys.stderr)
        return 2
    print(f"Report: {(args.out / 'report.html').resolve()}")
    print(f"Fits: {len(result['fits'])}; QC records: {len(result['qc'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
