#!/usr/bin/env python3

import argparse
import math
import shutil
import subprocess
import sys
from pathlib import Path


SKIP_EXIT_CODE = 77


def parse_args():
    parser = argparse.ArgumentParser(description="Run Viskores trace tests and optionally compare to MATLAB output.")
    parser.add_argument("--mode", choices=("generate", "compare"), required=True)
    parser.add_argument("--exe", required=True)
    parser.add_argument("--divertor", choices=("single", "circ"), required=True)
    parser.add_argument("--lines", required=True)
    parser.add_argument("--single-apar", required=True)
    parser.add_argument("--circ-apar", required=True)
    parser.add_argument("--reference-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--viskores-device", default="serial")
    parser.add_argument("--np-max", type=int, default=700)
    parser.add_argument("--tolerance", type=float, default=1.0e-8)
    return parser.parse_args()


def numeric_lines(lines_arg):
    out = []
    for token in lines_arg.split(","):
        token = token.strip()
        if token:
            out.append(int(token))
    if not out:
        raise ValueError("--lines did not contain any line numbers")
    return out


def line_number_from_output(path, divertor):
    prefix = f"ip_cxx.{divertor}."
    suffix = ".txt"
    name = path.name
    if not name.startswith(prefix) or not name.endswith(suffix):
        raise ValueError(f"Unexpected output filename: {name}")
    return int(name[len(prefix) : -len(suffix)])


def load_numeric_rows(path):
    rows = []
    for line in path.read_text().splitlines():
        values = []
        for token in line.split():
            try:
                values.append(float(token))
            except ValueError:
                pass
        if len(values) >= 5:
            rows.append(values)
    return rows


def check_inputs(args):
    exe = Path(args.exe)
    single_apar = Path(args.single_apar)
    circ_apar = Path(args.circ_apar)
    reference_dir = Path(args.reference_dir)

    missing = []
    if not exe.exists():
        missing.append(str(exe))
    if args.divertor == "single" and not single_apar.exists():
        missing.append(str(single_apar))
    if args.divertor == "circ" and not circ_apar.exists():
        missing.append(str(circ_apar))
    if args.mode == "compare" and not reference_dir.exists():
        missing.append(str(reference_dir))

    if missing:
        print("Skipping test because required files are missing:")
        for path in missing:
            print(f"  {path}")
        return False
    return True


def run_trace(args):
    output_dir = Path(args.output_dir)
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        args.exe,
        "--trace-engine",
        "viskores",
        "--viskores-device",
        args.viskores_device,
        "--viskores-output",
        "punctures",
        "--divertor",
        args.divertor,
        "--lines",
        args.lines,
        "--np-max",
        str(args.np_max),
        "--output-dir",
        str(output_dir),
    ]

    if args.divertor == "single":
        command += ["--single-apar", args.single_apar]
    else:
        command += ["--circ-apar", args.circ_apar]

    print("Running:")
    print(" ".join(command), flush=True)
    subprocess.run(command, check=True)


def combine_outputs(args):
    output_dir = Path(args.output_dir)
    files = [
        path
        for path in output_dir.glob(f"ip_cxx.{args.divertor}.*.txt")
        if ".TP." not in path.name and ".all." not in path.name
    ]
    files.sort(key=lambda path: line_number_from_output(path, args.divertor))

    combined = output_dir / f"ip_cxx.{args.divertor}.all.txt"
    with combined.open("w") as dst:
        wrote_header = False
        for path in files:
            lines = path.read_text().splitlines()
            if not lines:
                continue
            if not wrote_header:
                dst.write(lines[0] + "\n")
                wrote_header = True
            for line in lines[1:]:
                if line.strip():
                    dst.write(line + "\n")

    rows = max(0, sum(1 for _ in combined.open()) - 1)
    print(f"Combined puncture file: {combined}")
    print(f"Combined puncture rows: {rows}")


def compare_outputs(args):
    output_dir = Path(args.output_dir)
    reference_dir = Path(args.reference_dir)
    failures = 0

    print(f"\nComparing Viskores {args.divertor} output to MATLAB reference")
    print("line gen/ref rows maxAbsXYZ l2XYZ status")

    for line in numeric_lines(args.lines):
        generated = output_dir / f"ip_cxx.{args.divertor}.{line}.txt"
        reference = reference_dir / f"ip_matlab.{args.divertor}.{line}.txt"
        if not generated.exists() or not reference.exists():
            print(f"{line:4d} missing generated={generated.exists()} reference={reference.exists()} FAIL")
            failures += 1
            continue

        generated_rows = load_numeric_rows(generated)
        reference_rows = load_numeric_rows(reference)
        row_count = min(len(generated_rows), len(reference_rows))
        max_abs = 0.0
        sum_sq = 0.0
        count = 0

        for generated_row, reference_row in zip(generated_rows, reference_rows):
            for column in range(2, 5):
                diff = generated_row[column] - reference_row[column]
                max_abs = max(max_abs, abs(diff))
                sum_sq += diff * diff
                count += 1

        l2 = math.sqrt(sum_sq / count) if count else 0.0
        ok = len(generated_rows) == len(reference_rows) and max_abs < args.tolerance
        status = "PASS" if ok else "FAIL"
        if not ok:
            failures += 1

        print(f"{line:4d} {len(generated_rows):4d}/{len(reference_rows):4d} {max_abs:.12g} {l2:.12g} {status}")

    if failures:
        print(f"\nMATLAB comparison failed for {failures} case(s).")
        return 1

    print("\nMATLAB comparison passed.")
    return 0


def main():
    args = parse_args()
    if not check_inputs(args):
        return SKIP_EXIT_CODE

    run_trace(args)
    combine_outputs(args)

    if args.mode == "compare":
        return compare_outputs(args)

    return 0


if __name__ == "__main__":
    sys.exit(main())
