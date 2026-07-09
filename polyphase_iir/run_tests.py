import argparse
import json
import os
import random
import subprocess
import sys
from mpmath import mp
import numpy as np
from tabulate import tabulate

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from signal_mp import (
    butter_zpk_mp,
    ellip_zpk_mp,
    zpk2sos_mp,
    filter_polyphase_hybrid_mp,
    sosfilt_mp,
)
from design import (
    apply,
    design_polyphase_butterworth,
    design_polyphase_elliptic,
    generate_polyphase_cpp_struct,
    generate_polyphase_sos_cpp_struct,
    generate_polyphase_sos2_cpp_struct,
    format_cpp_array_2d,
)


def run_command(cmd, desc="command"):
    print(f"Running: {' '.join(cmd)}")
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"Error during {desc}!")
        print("STDOUT:", res.stdout)
        print("STDERR:", res.stderr)
        sys.exit(1)
    return res


def compute_reference_data(design_params, num_samples, out_path):
    filter_type = design_params[0]

    if filter_type == "butterworth":
        _, order, cutoff, M = design_params
        print(
            f"Designing Butterworth filter (Order={order}, Cutoff={cutoff}, M={M})..."
        )

        q_mp, a_sos_mp = design_polyphase_butterworth(
            order, cutoff, M, output="hybrid", as_float=False
        )
        polyphase_sos_mp = design_polyphase_butterworth(
            order, cutoff, M, output="sos", as_float=False
        )

        with mp.workdps(1000):
            serial_sos = zpk2sos_mp(*butter_zpk_mp(order, cutoff, btype="low"))

        ref_metadata = {
            "filter_type": "butterworth",
            "order": order,
            "cutoff": cutoff,
            "M": M,
        }

    elif filter_type == "elliptic":
        _, order, rp, rs, cutoff, M = design_params
        print(
            f"Designing Elliptic filter (Order={order}, Rp={rp}, Rs={rs}, Cutoff={cutoff}, M={M})..."
        )

        q_mp, a_sos_mp = design_polyphase_elliptic(
            order, rp, rs, cutoff, M, output="hybrid", as_float=False
        )
        polyphase_sos_mp = design_polyphase_elliptic(
            order, rp, rs, cutoff, M, output="sos", as_float=False
        )

        with mp.workdps(1000):
            serial_sos = zpk2sos_mp(*ellip_zpk_mp(order, rp, rs, cutoff, btype="low"))

        ref_metadata = {
            "filter_type": "elliptic",
            "order": order,
            "rp": rp,
            "rs": rs,
            "cutoff": cutoff,
            "M": M,
        }
    else:
        raise ValueError(f"Unsupported filter type: {filter_type}")

    with mp.workdps(100):
        random.seed(65429876)
        in_noise_mp = []
        for _ in range(M):
            in_noise_mp.append(
                [mp.mpf(random.uniform(-1.0, 1.0)) for _ in range(num_samples)]
            )

        print("Computing high-precision reference outputs...")
        out_noise_mp = filter_polyphase_hybrid_mp(q_mp, a_sos_mp, in_noise_mp)

        in_ir_mp = []
        out_ir_mp = []
        ir_length_per_branch = 65536 // M
        for i in range(M):
            sig = [mp.mpf(0) for _ in range(ir_length_per_branch)]
            sig[0] = mp.mpf(1)
            in_ir_mp.append(sig)

            sig = sosfilt_mp(polyphase_sos_mp[i], sig)
            out_ir_mp.append(sig)
        out_ir_mp = [v for w in zip(*out_ir_mp) for v in w]  # Interleave

        def _nstr(x, n=25):
            return mp.nstr(x, n)

        ref_data = {
            "serial_sos": apply(serial_sos, _nstr),
            "polyphase_sos": apply(polyphase_sos_mp, _nstr),
            "branches": apply(q_mp, _nstr),
            "denom": apply(a_sos_mp, _nstr),
            "inputs": apply(in_noise_mp, _nstr),
            "outputs": apply(out_noise_mp, _nstr),
            "impulse_response": apply(out_ir_mp, _nstr),
        }
        ref_data.update(ref_metadata)

        with open(out_path, "w") as f:
            json.dump(ref_data, f, indent=2)

    print(f"Reference data successfully generated at {out_path}\n")
    return ref_data


def generate_serial_sos_cpp_struct(sos, struct_name="ButterworthSos"):
    code = f"""template<typename T> struct {struct_name} {{
  static constexpr std::array<std::array<T, 5>, {len(sos)}> co{{{{
    {format_cpp_array_2d([s[:3] + s[4:6] for s in sos])}
  }}}};
}};"""
    return code


def generate_cpp_header(ref_data, header_path="coefficients.hpp"):
    with mp.workdps(100):
        code = "#pragma once\n#include <array>\n\nnamespace Coefficients {"

        code += "\n\n" + generate_polyphase_cpp_struct(
            apply(ref_data["branches"], lambda x: mp.mpf(x)),
            apply(ref_data["denom"], lambda x: mp.mpf(x)),
            "Hybrid",
        )

        code += "\n\n" + generate_polyphase_sos_cpp_struct(
            apply(ref_data["polyphase_sos"], lambda x: mp.mpf(x)),
            "PolyphaseSos",
        )

        code += "\n\n" + generate_polyphase_sos2_cpp_struct(
            apply(ref_data["polyphase_sos"], lambda x: mp.mpf(x)),
            "PolyphaseSos2",
        )

        code += "\n\n" + generate_serial_sos_cpp_struct(
            ref_data["serial_sos"], "SerialSos"
        )

        code += "\n\n}"

    with open(header_path, "w") as f:
        f.write(code)
    print(f"Code-generated coefficients header written to {header_path}")


def analyze_numerical_errors(ref_data, results_data, abs_tol):
    def measure(ref, test_vals, is_f32):
        test_outs = np.array([float(x) for x in test_vals], dtype=np.float64)

        # Absolute Error
        abs_err = np.abs(test_outs - ref)
        max_abs = np.max(abs_err)
        mean_abs = np.mean(abs_err)

        # Relative Error without Tol filtering (can skew heavily when outputs cross 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_err_no_tol = abs_err / np.abs(ref)
            rel_err_no_tol_clean = rel_err_no_tol[np.isfinite(rel_err_no_tol)]
        max_rel_no_tol = (
            np.max(rel_err_no_tol_clean) if len(rel_err_no_tol_clean) > 0 else 0.0
        )
        mean_rel_no_tol = (
            np.mean(rel_err_no_tol_clean) if len(rel_err_no_tol_clean) > 0 else 0.0
        )

        # Relative Error with absolute tolerance filter
        valid_mask = np.abs(ref) > abs_tol
        if np.any(valid_mask):
            rel_err_with_tol = abs_err[valid_mask] / np.abs(ref[valid_mask])
            max_rel_with_tol = np.max(rel_err_with_tol)
            mean_rel_with_tol = np.mean(rel_err_with_tol)
        else:
            max_rel_with_tol = 0.0
            mean_rel_with_tol = 0.0

        # ULP Error
        if is_f32:
            ref_f32 = np.array(ref, dtype=np.float32)
            ulp_sizes = np.abs(np.spacing(ref_f32)).astype(np.float64)
        else:
            ulp_sizes = np.abs(np.spacing(ref))

        # Avoid zero division inside subnormal boundaries
        ulp_sizes[ulp_sizes == 0] = (
            np.finfo(np.float32).tiny if is_f32 else np.finfo(np.float64).tiny
        )

        ulp_err = abs_err / ulp_sizes
        max_ulp = np.max(ulp_err)
        mean_ulp = np.mean(ulp_err)

        return {
            "abs_max": max_abs,
            "abs_mean": mean_abs,
            "rel_no_tol_max": max_rel_no_tol,
            "rel_no_tol_mean": mean_rel_no_tol,
            "rel_with_tol_max": max_rel_with_tol,
            "rel_with_tol_mean": mean_rel_with_tol,
            "ulp_max": max_ulp,
            "ulp_mean": mean_ulp,
        }

    # Render Report Table
    print("")
    filter_type = ref_data.get("filter_type", "butterworth")
    if filter_type == "elliptic":
        title_str = (
            f"## Filter={filter_type}, Order={ref_data['order']}, "
            f"Rp={ref_data.get('rp')}, Rs={ref_data.get('rs')}, "
            f"Cutoff={ref_data['cutoff']}, M={ref_data['M']}"
        )
    else:
        title_str = (
            f"## Filter={filter_type}, Order={ref_data['order']}, "
            f"Cutoff={ref_data['cutoff']}, M={ref_data['M']}"
        )
    print(title_str)
    print("Results are shown as Max / Mean.")
    print("")

    ref_noise = np.array([float(x) for x in ref_data["outputs"]], dtype=np.float64)
    ref_ir = np.array(
        [float(x) for x in ref_data["impulse_response"]], dtype=np.float64
    )

    f64_rows = []
    f32_rows = []
    for name, output_list in results_data.items():
        is_f32 = "float" in name or "f32" in name

        impl_title = name.replace("_", " ").title()

        if "_noise" in name:
            errs = measure(ref_noise, output_list, is_f32)
        else:
            if ref_ir is None:
                ref_ir = np.zeros(len(output_list))
                ref_ir[0] = 1
            errs = measure(ref_ir, output_list, is_f32)

        abs_str = f"{errs['abs_max']:.4e} / {errs['abs_mean']:.4e}"
        no_tol_str = f"{errs['rel_no_tol_max']:.4e} / {errs['rel_no_tol_mean']:.4e}"
        with_tol_str = (
            f"{errs['rel_with_tol_max']:.4e} / {errs['rel_with_tol_mean']:.4e}"
        )
        ulp_str = f"{errs['ulp_max']:>.2e} / {errs['ulp_mean']:>.2e}"

        row = [impl_title, abs_str, no_tol_str, with_tol_str, ulp_str]
        if is_f32:
            f32_rows.append(row)
        else:
            f64_rows.append(row)

    headers = [
        "Implementation",
        "Abs Error",
        "Rel Error, Tol=0",
        f"Rel Error, Tol={abs_tol:.2e}",
        "ULP",
    ]

    table_data = f64_rows + f32_rows
    table_str = tabulate(table_data, headers=headers, tablefmt="pipe")
    print(table_str)
    print("")


def print_benchmark_results(benchmarks):
    print("\n## Benchmarks\n")
    headers = [
        "Implementation",
        "Time (ms)",
        "ns / Sample",
        "Throughput (MSps)",
        "nSample",
    ]

    f64_rows = []
    f32_rows = []
    for name, stats in benchmarks.items():
        is_f32 = "float" in name or "f32" in name
        impl_title = name.replace("_", " ").title()
        time_ms = stats.get("time_ms", 0.0)
        ns_per_sample = stats.get("ns_per_sample", 0.0)
        num_samples = stats.get("num_samples", 1)
        throughput_msps = (num_samples / time_ms) / 1000.0 if time_ms > 0 else 0.0

        row = [impl_title, time_ms, ns_per_sample, throughput_msps, num_samples]
        if is_f32:
            f32_rows.append(row)
        else:
            f64_rows.append(row)

    table_data = f64_rows + f32_rows
    # Using 'floatfmt' to apply identical decimal formatting to all float columns
    table_str = tabulate(table_data, headers=headers, tablefmt="pipe", floatfmt=".2f")
    print(table_str)


def main():
    parser = argparse.ArgumentParser(
        description="Test framework for PolyphaseIir numerical errors."
    )
    parser.add_argument(
        "--design",
        type=str,
        help=(
            "Design parameters.\n"
            "  'butterworth,order,cutoff,M' (e.g., 'butterworth,8,0.2,4')\n"
            "  'elliptic,order,rp,rs,cutoff,M' (e.g., 'elliptic,8,0.1,60,0.2,4')"
        ),
    )
    parser.add_argument(
        "--ref",
        type=str,
        default="reference.json",
        help="Path to save or read the reference data JSON",
    )
    parser.add_argument(
        "--results",
        type=str,
        default="results.json",
        help="Path to save the results data JSON",
    )
    parser.add_argument(
        "--samples",
        type=int,
        default=48000,
        help="Number of signal samples to generate",
    )
    parser.add_argument(
        "--compiler", type=str, default="g++", help="C++ compiler executable"
    )
    parser.add_argument(
        "--abs-tol",
        type=float,
        default=np.finfo(np.float32).eps,
        help="Absolute tolerance threshold for relative error filtering (default: f64 machine epsilon)",
    )
    parser.add_argument(
        "--benchmark",
        action="store_true",
        help="Enable and run execution speed benchmarks",
    )

    args = parser.parse_args()

    if os.path.exists(args.ref):
        print(f"Loading existing reference file from {args.ref}")
        with open(args.ref) as f:
            ref_data = json.load(f)
    elif args.design:
        try:
            parts = args.design.split(",")
            filter_type = parts[0].strip().lower()

            if filter_type in ("butterworth", "butter"):
                order = int(parts[1])
                cutoff = float(parts[2])
                M = int(parts[3])
                design_params = ("butterworth", order, cutoff, M)
            elif filter_type in ("elliptic", "ellip"):
                order = int(parts[1])
                rp = float(parts[2])
                rs = float(parts[3])
                cutoff = float(parts[4])
                M = int(parts[5])
                design_params = ("elliptic", order, rp, rs, cutoff, M)
            else:
                raise ValueError(f"Unknown filter type: '{filter_type}'")
        except Exception as e:
            print(f"Error parsing --design parameters: {e}")
            print("Use format:")
            print("  butterworth,order,cutoff,M       (e.g., 'butterworth,8,0.2,4')")
            print(
                "  elliptic,order,rp,rs,cutoff,M    (e.g., 'elliptic,8,0.1,60,0.2,4')"
            )
            sys.exit(1)

        ref_data = compute_reference_data(design_params, args.samples, args.ref)
    else:
        print(
            f"No reference data found at {args.ref}. Designing default filter 'butterworth,8,0.2,4'..."
        )
        default_design = ("butterworth", 8, 0.2, 4)
        ref_data = compute_reference_data(default_design, args.samples, args.ref)

    generate_cpp_header(ref_data)

    is_cl = args.compiler.lower() == "cl" or args.compiler.lower() == "cl.exe"
    if is_cl:
        compile_cmd = [
            args.compiler,
            "/std:c++20",
            "/O2",
            "/EHsc",
            "/arch:AVX2",
            "test_runner.cpp",
            "/Fetest_runner",
        ]
        if args.benchmark:
            compile_cmd.append("/DENABLE_BENCHMARK=1")
    else:
        compile_cmd = [
            args.compiler,
            "-std=c++20",
            "-O3",
            "-mfma",
            "-march=native",
            "test_runner.cpp",
            "-o",
            "test_runner",
        ]
        if args.benchmark:
            compile_cmd.append("-DENABLE_BENCHMARK=1")

    run_command(compile_cmd, "C++ test runner compilation")
    run_command(
        [f"./test_runner{'.exe' if is_cl else ''}", args.ref, args.results],
        "C++ test runner execution",
    )

    with open(args.results) as f:
        results_data = json.load(f)

    if "benchmarks" in results_data:
        print_benchmark_results(results_data["benchmarks"])
        results_data = {k: v for k, v in results_data.items() if k != "benchmarks"}

    analyze_numerical_errors(ref_data, results_data, args.abs_tol)


if __name__ == "__main__":
    main()
