"""
plot_results.py

Visualization utility for Polyphase IIR filter test results.
- Extracts the frequency responses directly from the implemented impulse responses
  (untructated) to accurately capture magnitude, phase, and group delay.
- Computes the numerical error comparison for the quality scatter plot using
  the filtered white noise outputs.
- Compares magnitude, phase, group delay, and impulse response against the
  mathematical reference curves.
- Plots filter quality as a D scatter plot: Numerical Error vs. Execution Speed.
- Includes intuitive configuration switches to easily toggle specific filters.

Usage:
  1. Generate the test results (e.g., using a benchmark run):
     python run_tests.py --compiler g++ --design 8,0.03125,4 --samples 1000 --benchmark
  2. Run the visualization script:
     python plot_results.py
"""

import argparse
import json
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

# =====================================================================
# FILTER PLOTTING CONFIGURATION (Toggles)
# =====================================================================
# Toggle global precision classes
PLOT_F64 = True
PLOT_F32 = True

# Individual filter implementation switches
PLOT_CONFIG = {
    "reference": True,  # Plots the exact mathematical Butterworth reference
    "tdf2_f64": True,
    "tdf2_f32": True,
    "hybrid_simple_f64": True,
    "hybrid_simple_f32": True,
    "hybrid_kahan_f64": True,
    "hybrid_kahan_f32": True,
    "sos_simple_f64": True,
    "sos_simple_f32": True,
    "sos_kahan_f64": True,
    "sos_kahan_f32": True,
}

# =====================================================================
# VISUAL STYLING
# =====================================================================
FILTER_STYLES = {
    "tdf2_f64": {
        "color": "#1f77b4",
        "marker": "o",
        "label": "TDF2 (f64)",
        "linestyle": "-",
    },
    "hybrid_simple_f64": {
        "color": "#aec7e8",
        "marker": "s",
        "label": "Hybrid Simple (f64)",
        "linestyle": "-",
    },
    "hybrid_kahan_f64": {
        "color": "#ff7f0e",
        "marker": "D",
        "label": "Hybrid Kahan (f64)",
        "linestyle": "-",
    },
    "sos_simple_f64": {
        "color": "#ffbb78",
        "marker": "^",
        "label": "SOS Simple (f64)",
        "linestyle": "-",
    },
    "sos_kahan_f64": {
        "color": "#2ca02c",
        "marker": "v",
        "label": "SOS Kahan (f64)",
        "linestyle": "-",
    },
    "tdf2_f32": {
        "color": "#1f77b4",
        "marker": "o",
        "label": "TDF2 (f32)",
        "linestyle": "--",
        "fillstyle": "none",
    },
    "hybrid_simple_f32": {
        "color": "#aec7e8",
        "marker": "s",
        "label": "Hybrid Simple (f32)",
        "linestyle": "--",
        "fillstyle": "none",
    },
    "hybrid_kahan_f32": {
        "color": "#ff7f0e",
        "marker": "D",
        "label": "Hybrid Kahan (f32)",
        "linestyle": "--",
        "fillstyle": "none",
    },
    "sos_simple_f32": {
        "color": "#ffbb78",
        "marker": "^",
        "label": "SOS Simple (f32)",
        "linestyle": "--",
        "fillstyle": "none",
    },
    "sos_kahan_f32": {
        "color": "#2ca02c",
        "marker": "v",
        "label": "SOS Kahan (f32)",
        "linestyle": "--",
        "fillstyle": "none",
    },
}

# Set style if available
try:
    plt.style.use("seaborn-v0_8-whitegrid")
except Exception:
    try:
        plt.style.use("ggplot")
    except Exception:
        pass


def relative_error(ref, approx, dtype):
    ref_arr = np.asarray(ref, dtype=dtype)
    approx_arr = np.asarray(approx, dtype=dtype)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = np.abs((ref_arr - approx_arr) / ref_arr, dtype=dtype)
    return np.nan_to_num(error, nan=0.0)


def main():
    parser = argparse.ArgumentParser(
        description="Visualize Polyphase IIR filter test results."
    )
    parser.add_argument(
        "--ref", type=str, default="reference.json", help="Path to reference.json"
    )
    parser.add_argument(
        "--results", type=str, default="results.json", help="Path to results.json"
    )
    parser.add_argument(
        "--save",
        type=str,
        default="",
        help="Save plots as images instead of displaying them",
    )
    args = parser.parse_args()

    # Load data
    if not os.path.exists(args.ref):
        print(
            f"Error: Reference file '{args.ref}' not found. Please run run_tests.py first."
        )
        sys.exit(1)
    if not os.path.exists(args.results):
        print(
            f"Error: Results file '{args.results}' not found. Please run run_tests.py first."
        )
        sys.exit(1)

    with open(args.ref) as f:
        ref_data = json.load(f)
    with open(args.results) as f:
        results_data = json.load(f)

    M = int(ref_data["M"])
    cutoff = float(ref_data["cutoff"])
    order = int(ref_data["order"])
    benchmarks = results_data.get("benchmarks", {})

    worN = np.linspace(0, min(1.1 * np.pi / M, np.pi), 2048)

    ref_ir = np.array(
        [float(x) for x in ref_data["impulse_response"]], dtype=np.float64
    )
    w_ref, h_ref = signal.freqz(ref_ir, 1, worN=worN)
    gain_ref = 20 * np.log10(np.abs(h_ref) + 1e-15)
    phase_ref = np.unwrap(np.angle(h_ref))
    _, gd_ref = signal.group_delay((ref_ir, 1), w=worN)

    ir_length = min(len(ref_ir), 2048)

    # Dictionary to keep extracted results
    extracted = {}
    errors = {}

    ref_outputs = np.array([float(x) for x in ref_data["outputs"]], dtype=np.float64)

    # Process filters based on configured base implementation names
    base_names = [
        "tdf2_f64",
        "hybrid_simple_f64",
        "hybrid_kahan_f64",
        "sos_simple_f64",
        "sos_kahan_f64",
        "tdf2_f32",
        "hybrid_simple_f32",
        "hybrid_kahan_f32",
        "sos_simple_f32",
        "sos_kahan_f32",
    ]

    for name in base_names:
        # Apply group and individual configurations
        if "f32" in name and not PLOT_F32:
            continue
        if "f64" in name and not PLOT_F64:
            continue
        if name in PLOT_CONFIG and not PLOT_CONFIG[name]:
            continue

        # 1. Retrieve the actual (untruncated) impulse response for frequency response plotting
        ir_key = f"{name}_ir"
        if ir_key in results_data:
            h_est = np.array([float(x) for x in results_data[ir_key]], dtype=np.float64)
            is_diverged = not np.all(np.isfinite(h_est)) or np.max(np.abs(h_est)) > 10.0
            extracted[name] = (h_est, is_diverged)
        else:
            continue

        # 2. Retrieve noise outputs for quality scatter plot error calculations
        noise_key = f"{name}_noise"
        if noise_key in results_data:
            dtype = np.float32 if "f32" in name else np.float64
            y_test = np.array(
                [float(val) for val in results_data[noise_key]], dtype=dtype
            )
            if len(y_test) == len(ref_outputs):
                abs_err = np.abs(y_test - ref_outputs)
                errors[name] = np.max(abs_err)
            else:
                errors[name] = np.nan
        else:
            errors[name] = np.nan

    # =====================================================================
    # FIGURE 1: FREQUENCY RESPONSES
    # =====================================================================
    fig1, axs = plt.subplots(2, 2, figsize=(15, 8))
    fig1.suptitle(
        f"Frequency Response Analysis\n"
        f"(Order={order}, Cutoff={cutoff}, Decimation M={M})",
        fontsize=14,
        weight="bold",
    )

    ax_gain, ax_phase = axs[0, 0], axs[0, 1]
    ax_gd, ax_ir = axs[1, 0], axs[1, 1]

    # Plot Reference Curves
    if PLOT_CONFIG.get("reference", True):
        ax_gain.plot(
            w_ref / np.pi,
            gain_ref,
            label="Reference (mpmath)",
            color="black",
            linestyle="-",
            linewidth=2.0,
            zorder=2,
        )
        ax_phase.plot(
            w_ref / np.pi,
            phase_ref,
            label="Reference (mpmath)",
            color="black",
            linestyle="-",
            linewidth=2.0,
            zorder=2,
        )
        ax_gd.plot(
            w_ref / np.pi,
            gd_ref,
            label="Reference (mpmath)",
            color="black",
            linestyle="-",
            linewidth=2.0,
            zorder=2,
        )
        ax_ir.plot(
            np.arange(min(ir_length, len(ref_ir))),
            ref_ir[:ir_length],
            label="Reference (mpmath)",
            color="black",
            linestyle="-",
            linewidth=2.0,
            zorder=2,
        )

    # Plot each enabled filter
    for name, (h_est, is_diverged) in extracted.items():
        style = FILTER_STYLES.get(name, {})
        label = style.get("label", name)
        color = style.get("color", None)
        ls = style.get("linestyle", "-")

        if is_diverged:
            print(f"Warning: Filter '{name}' shows instability/divergence.")
            ax_ir.plot(
                np.arange(min(ir_length, len(h_est))),
                h_est[:ir_length],
                label=f"{label} (DIVERGED)",
                color=color,
                linestyle=ls,
            )
            continue

        # Compute frequency parameters from the full untruncated impulse response
        w, h = signal.freqz(h_est, 1, worN=worN)
        gain = 20 * np.log10(np.abs(h) + 1e-15)
        phase = np.unwrap(np.angle(h))
        _, gd = signal.group_delay((h_est, 1), w=worN)

        # Plot curves
        ax_gain.plot(w / np.pi, gain, label=label, color=color, linestyle=ls)
        ax_phase.plot(w / np.pi, phase, label=label, color=color, linestyle=ls)
        ax_gd.plot(w / np.pi, gd, label=label, color=color, linestyle=ls)
        ax_ir.plot(
            np.arange(min(ir_length, len(h_est))),
            h_est[:ir_length],
            label=label,
            color=color,
            linestyle=ls,
        )

    # Format Subplots (Limits, Reference Lines, Labels)
    for ax in [ax_gain, ax_phase, ax_gd]:
        ax.axvline(cutoff, color="red", linestyle=":", alpha=0.8, label="Filter Cutoff")
        ax.axvline(
            1.0 / M,
            color="purple",
            linestyle="--",
            alpha=0.8,
            label="Decimation Nyquist (1/M)",
        )
        ax.set_xlabel("Normalized Frequency (rad/sample)")

    ax_gain.set_ylabel("Magnitude (dB)")
    ax_gain.set_title("Magnitude Response")
    # ax_gain.set_ylim([-90, 5])

    ax_phase.set_ylabel("Phase (radians)")
    ax_phase.set_title("Unwrapped Phase Response")

    ax_gd.set_ylabel("Group Delay (samples)")
    ax_gd.set_title("Group Delay Response")

    # Bind group delay y-axis to reference range to ignore numerical spikes near stopband
    valid_gd_ref = gd_ref[gain_ref > -80.0]
    if len(valid_gd_ref) > 0:
        gd_ref_max = np.nanmax(valid_gd_ref)
        gd_ref_min = np.nanmin(valid_gd_ref)
        margin = max(1.0, 0.2 * (gd_ref_max - gd_ref_min))
        ax_gd.set_ylim([max(0, gd_ref_min - margin), gd_ref_max + margin])

    ax_ir.set_xlabel("Time Index (samples)")
    ax_ir.set_ylabel("Amplitude")
    ax_ir.set_title("Impulse Response (Trimmed)")
    ax_ir.axhline(0, color="gray", linestyle=":", alpha=0.5)

    # Scale IR view window based on stable filters
    stable_h_vals = [h[:ir_length] for h, div in extracted.values() if not div]
    if stable_h_vals:
        max_stable = np.max(np.abs(stable_h_vals))
        ax_ir.set_ylim([-1.2 * max_stable, 1.2 * max_stable])

    for ax in axs.flat:
        ax.legend(fontsize=8, loc="best")
        ax.grid(True, which="both", linestyle=":", alpha=0.5)

    plt.tight_layout()

    # =====================================================================
    # FIGURE 2: QUALITY SCATTER PLOTS (Accuracy vs. Speed - Split f32 / f64)
    # =====================================================================
    fig2, (ax_f32, ax_f64) = plt.subplots(1, 2, figsize=(12, 4))
    fig2.suptitle(
        "Filter Implementation Quality (Bottom-left is better)",
        fontsize=14,
        weight="bold",
    )

    has_benchmarks = len(benchmarks) > 0

    # Define unique configuration settings for each subplot
    axes_config = {
        "f32": {
            "ax": ax_f32,
            "title": "f32 Single Precision Implementation",
            "epsilons": [
                (np.finfo(np.float32).eps, "f32 Epsilon (~1.19e-07)", "green")
            ],
        },
        "f64": {
            "ax": ax_f64,
            "title": "f64 Double Precision Implementation",
            "epsilons": [
                (np.finfo(np.float64).eps, "f64 Epsilon (~2.22e-16)", "blue"),
            ],
        },
    }

    # Set up basic structure, titles, axes labels, and background helpers
    for key, cfg in axes_config.items():
        ax = cfg["ax"]
        ax.set_title(cfg["title"], fontsize=11, weight="bold")
        ax.set_yscale("log")
        ax.set_ylabel("Numerical Precision: Max Absolute Error")

        if has_benchmarks:
            ax.set_xlabel("Execution Time per Sample (ns / sample)")
        else:
            ax.set_xlabel("Execution Time (Benchmark Data Not Available)")
            ax.set_xticks([])
            ax.text(
                0.5,
                0.5,
                "Notice: Speed metrics are unavailable.\n"
                "To capture execution speed, run the tests using:\n"
                "python run_tests.py --benchmark ...",
                color="red",
                fontsize=10,
                weight="bold",
                ha="center",
                va="center",
                transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="red", alpha=0.8),
            )

        ax.grid(True, which="both", linestyle=":", alpha=0.5)

    # Plot metrics to their respective target axes
    for name in errors:
        style = FILTER_STYLES.get(name, {})
        label = style.get("label", name)
        color = style.get("color", "#7f7f7f")
        marker = style.get("marker", "o")
        fillstyle = style.get("fillstyle", "full")

        max_err = errors[name]
        is_diverged = extracted[name][1] if name in extracted else True

        # Handle Y coordinate (Accuracy / Error)
        if is_diverged or not np.isfinite(max_err):
            plot_err = 10.0  # Force to top margin to indicate divergence
            text_color = "red"
        else:
            plot_err = max_err
            text_color = "black"

        # Handle X coordinate (Execution Time per Sample)
        if has_benchmarks and name in benchmarks:
            plot_time = benchmarks[name]["ns_per_sample"]
        else:
            plot_time = 1.0  # Dummy value

        facecolor = "none" if fillstyle == "none" else color

        # Select the target subplot axis based on the data precision type
        if "f32" in name:
            ax = ax_f32
        elif "f64" in name:
            ax = ax_f64
        else:
            continue

        ax.scatter(
            plot_time,
            plot_err,
            facecolors=facecolor,
            edgecolors=color,
            marker=marker,
            s=120,
            label=label,
            zorder=3,
        )

    # Add reference machine epsilons and adjust vertical boundaries dynamically
    for key, cfg in axes_config.items():
        ax = cfg["ax"]

        # Draw epsilon guides
        for eps_val, eps_label, eps_color in cfg["epsilons"]:
            ax.axhline(
                eps_val, color=eps_color, linestyle=":", alpha=0.5, label=eps_label
            )

        key_errors = [
            errors[name] for name in errors if key in name and np.isfinite(errors[name])
        ]
        if key_errors:
            min_err = min(key_errors)
            max_err = max(key_errors)
            ax.set_ylim([min_err * 0.9, max_err * 1.1])

        ax.legend()

    plt.tight_layout()

    if len(args.save) > 0:
        os.makedirs("img", exist_ok=True)
        fig1.savefig(f"img/response_{args.save}.svg")
        fig2.savefig(f"img/quality_{args.save}.svg")
    else:
        plt.show()


if __name__ == "__main__":
    main()
