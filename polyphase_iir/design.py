import numpy as np
import matplotlib.pyplot as plt
import copy
from mpmath import mp
from scipy import signal
from signal_mp import butter_zpk_mp, ellip_zpk_mp, tf2sos_mp


def apply(data, fn, deepcopy=True):
    """
    Apply `fn` recusrively to the nested lists.
    """
    if deepcopy:
        res = copy.deepcopy(data)
    else:
        res = data
    stack = [res]
    while stack:
        current = stack.pop()
        if isinstance(current, list):
            for i in range(len(current)):
                if isinstance(current[i], list):
                    stack.append(current[i])
                else:
                    current[i] = fn(current[i])
    return res


def design_polyphase_iir(
    z_mp,
    p_mp,
    k_mp,
    M: int,
    output: str = "ba",
    as_float: bool = False,
    workdps: int | None = None,
):
    """
    Decomposes an arbitrary IIR filter (given by its zeros, poles, and gain)
    for polyphase down-sampling.

    Parameters:
    -----------
    z_mp : list or iterable
        Zeros of the filter as mpmath numbers.
    p_mp : list or iterable
        Poles of the filter as mpmath numbers.
    k_mp : mpmath number
        System gain.
    M : int
        Down-sampling factor.
    output : str, optional
        Format of the output filter: "ba", "sos", "sos2", or "hybrid".
    as_float : bool, optional
        Cast the coefficients to Python float.
    workdps : int, optional
        Working precision in decimal places (dps). If None, the current
        ambient mpmath precision context is used.

    Returns:
    --------
    tuple or list
        The returned structure depends on the value of the `output` parameter:
        * If "ba" (default):
          Returns a tuple `(q_polyphase, a_low_poly)` where:
            - `q_polyphase` (list of lists): Numerator coefficients for each of the `M` polyphase branches.
            - `a_low_poly` (list): Denominator coefficients shared by all branches (expressed in powers of z^-M).
        * If "sos":
          Returns `sos_polyphase` (list of lists of lists): A list of length `M` containing the Second-Order Section (SOS) representations for each polyphase branch.
        * If "sos2":
          Returns a tuple `(b_all, a_all)` optimized for shared denominators:
            - `b_all` (list of lists of lists): For each of the `M` branches, a list
              of numerator coefficients `[b0, b1, b2]` for each SOS section.
            - `a_all` (list of lists): A single list of denominator coefficients
              `[a1, a2]` for each SOS section (normalized $a0=1$), shared by all branches.
        * If "hybrid":
          Returns a tuple `(q_polyphase, sos_sections)` where:
            - `q_polyphase` (list of lists): Numerator coefficients for each of the `M` polyphase branches.
            - `sos_sections` (list of lists): Denominator coefficients `[a1, a2]` for each pole pair.

        Coefficients are returned as standard Python floats if `as_float` is True; otherwise, they are `mpmath` types.
    """
    if output not in ("ba", "sos", "sos2", "hybrid"):
        raise ValueError("output must be one of 'ba', 'sos', 'sos2', or 'hybrid'")

    with mp.workdps(workdps):

        def poly_from_roots(roots):
            """
            Expands (x - r1)(x - r2)... into polynomial coefficients.
            Returns coefficients in descending order of powers (highest power first).
            """
            c = [mp.mpf(1)]
            for r in roots:
                c += [0]
                for i in range(len(c) - 1, 0, -1):
                    c[i] -= r * c[i - 1]
            return c

        def convolve(a, b):
            len_a = len(a)
            len_b = len(b)
            out = [mp.mpf(0.0)] * (len_a + len_b - 1)
            for i in range(len_a):
                for j in range(len_b):
                    out[i + j] += a[i] * b[j]
            return out

        b_poly = [mp.re(x * k_mp) for x in poly_from_roots(z_mp)]

        p_new = [p**M for p in p_mp]
        a_low_poly = [mp.re(x) for x in poly_from_roots(p_new)]

        s_poly = [mp.mpc(1.0)]
        for pi in p_mp:
            S_i = [pi**k for k in range(M)]
            s_poly = convolve(s_poly, S_i)
        s_poly = [mp.re(x) for x in s_poly]

        q_poly = [mp.re(x) for x in convolve(b_poly, s_poly)]

        q_polyphase = []
        for k in range(M):
            q_polyphase.append(q_poly[k::M])

        if output == "ba":
            if as_float:
                q_polyphase = apply(q_polyphase, float)
                a_low_poly = apply(a_low_poly, float)
            return q_polyphase, a_low_poly

        elif output == "sos":
            sos_polyphase = []
            for k in range(M):
                q_k = q_poly[k::M]
                sos_section = tf2sos_mp(q_k, a_low_poly)
                sos_polyphase.append(sos_section)

            if as_float:
                sos_polyphase = apply(sos_polyphase, float)
            return sos_polyphase

        elif output == "sos2":
            sos_polyphase = []
            for k in range(M):
                q_k = q_poly[k::M]
                sos_section = tf2sos_mp(q_k, a_low_poly)
                sos_polyphase.append(sos_section)

            b_all = []
            a_all = []
            nSections = len(sos_polyphase[0]) if len(sos_polyphase) > 0 else 0

            for phase in sos_polyphase:
                b_phase = []
                for sec in phase:
                    b0, b1, b2, a0, a1, a2 = sec
                    b_phase.append([b0 / a0, b1 / a0, b2 / a0])
                b_all.append(b_phase)

            if nSections > 0:
                for s_idx in range(nSections):
                    b0, b1, b2, a0, a1, a2 = sos_polyphase[0][s_idx]
                    a_all.append([a1 / a0, a2 / a0])

            if as_float:
                b_all = apply(b_all, float)
                a_all = apply(a_all, float)
            return b_all, a_all

        elif output == "hybrid":
            sos_sections = []
            pool = list(p_new)
            while len(pool) > 0:
                p1 = pool.pop(0)

                if abs(mp.im(p1)) < 1e-20:
                    a1 = -mp.re(p1)
                    a2 = mp.mpf(0.0)
                    sos_sections.append([a1, a2])
                    continue

                best_idx = -1
                best_err = mp.mpf("inf")

                for i, p2 in enumerate(pool):
                    err = abs(p1 - mp.conj(p2))
                    if err < best_err:
                        best_err = err
                        best_idx = i

                if best_idx != -1:
                    p2 = pool.pop(best_idx)
                    a1 = -2 * mp.re(p1)
                    a2 = abs(p1) ** 2
                    sos_sections.append([a1, a2])
                else:
                    a1 = -mp.re(p1)
                    a2 = mp.mpf(0.0)
                    sos_sections.append([a1, a2])

            if as_float:
                q_polyphase = apply(q_polyphase, float)
                sos_sections = apply(sos_sections, float)
            return q_polyphase, sos_sections


def design_polyphase_butterworth(
    order: int,
    cutoff: float,
    M: int,
    output: str = "ba",
    as_float: bool = False,
    workdps: int | None = 2000,
):
    """
    Designs a Butterworth filter and decomposes it for polyphase down-sampling.

    Parameters:
    -----------
    order : int
        Order of the Butterworth filter.
    cutoff : float
        Normalized cutoff frequency (0 < cutoff < 1).
    M : int
        Down-sampling factor.
    output : str, optional
        Format of the output filter. Can be "ba", "sos", or "hybrid".
    as_float : bool, optional
        Cast the coefficients to Python float.
    workdps : int, optional
        Working precision in dps. If None, it is estimated.

    Returns:
    --------
    tuple or list
        The structured output depending on the `output` parameter (refer to `design_polyphase_iir` for details):
        * If "ba": A tuple `(q_polyphase, a_low_poly)`.
        * If "sos": A list `sos_polyphase` containing SOS representations for each branch.
        * If "sos2": A tuple `(b_all, a_all)` separating numerator and shared denominator sections.
        * If "hybrid": A tuple `(q_polyphase, sos_sections)`.

        Coefficients are returned as standard Python floats if `as_float` is True; otherwise, they are high-precision `mpmath` types.
    """
    if output not in ("ba", "sos", "sos2", "hybrid"):
        raise ValueError("output must be one of 'ba', 'sos', 'sos2', or 'hybrid'")

    if workdps is None:
        f_min = min(cutoff, 1.0 - cutoff)
        log_fc = np.log10(f_min)

        slope = -0.8 * log_fc - 0.2
        intercept = 0.3 * log_fc - 0.5
        est_peak_1a = slope * order + intercept

        upper_bound = order * np.log10(M)

        est_log_max_s = max(0.0, min(upper_bound, est_peak_1a))
        workdps = int(np.ceil(est_log_max_s)) + 25
        workdps = max(25, workdps)

    with mp.workdps(workdps):
        z_mp, p_mp, k_mp = butter_zpk_mp(order, cutoff, btype="low")
    return design_polyphase_iir(
        z_mp, p_mp, k_mp, M, output=output, as_float=as_float, workdps=workdps
    )


def design_polyphase_elliptic(
    order: int,
    rp: float,
    rs: float,
    cutoff: float,
    M: int,
    output: str = "ba",
    as_float: bool = False,
    workdps: int | None = 2000,
):
    """
    Designs an Elliptic filter and decomposes it for polyphase down-sampling.

    Parameters:
    -----------
    order : int
        Order of the Elliptic filter.
    rp : float
        Passband ripple in dB.
    rs : float
        Stopband attenuation in dB.
    cutoff : float
        Normalized cutoff frequency (0 < cutoff < 1).
    M : int
        Down-sampling factor.
    output : str, optional
        Format of the output filter. Can be "ba", "sos", or "hybrid".
    as_float : bool, optional
        Cast the coefficients to Python float.
    workdps : int, optional
        Working precision in dps. If None, it is estimated.

    Returns:
    --------
    tuple or list
        The structured output depending on the `output` parameter (refer to `design_polyphase_iir` for details):
        * If "ba": A tuple `(q_polyphase, a_low_poly)`.
        * If "sos": A list `sos_polyphase` containing SOS representations for each branch.
        * If "sos2": A tuple `(b_all, a_all)` separating numerator and shared denominator sections.
        * If "hybrid": A tuple `(q_polyphase, sos_sections)`.

        Coefficients are returned as standard Python floats if `as_float` is True; otherwise, they are high-precision `mpmath` types.
    """
    if output not in ("ba", "sos", "sos2", "hybrid"):
        raise ValueError("output must be one of 'ba', 'sos', 'sos2', or 'hybrid'")

    if workdps is None:
        f_min = min(cutoff, 1.0 - cutoff)
        log_fc = np.log10(f_min)

        slope = -0.8 * log_fc - 0.2
        intercept = 0.3 * log_fc - 0.5
        est_peak_1a = slope * order + intercept

        upper_bound = order * np.log10(M)

        est_log_max_s = max(0.0, min(upper_bound, est_peak_1a))
        workdps = int(np.ceil(est_log_max_s)) + 25
        workdps = max(25, workdps)

    with mp.workdps(workdps):
        z_mp, p_mp, k_mp = ellip_zpk_mp(order, rp, rs, cutoff, btype="low")
    return design_polyphase_iir(
        z_mp, p_mp, k_mp, M, output=output, as_float=as_float, workdps=workdps
    )


def check_stability(order, M, precision="float64"):
    """
    Checks if the Polyphase Butterworth design will be numerically stable.
    Formula derived from regression of coefficient explosion.
    """

    if order <= 8:
        return True

    est_log_mag = 0.015 * M * ((order - 8) ** 2) - 10

    limit = 14.0 if precision == "float64" else 5.0

    if est_log_mag > limit:
        print(f"WARNING: Numerical Instability Predicted!")
        print(f"  Params: Order={order}, Decimation={M}")
        print(f"  Est. Coefficient Range: +/- 1e{est_log_mag:.1f}")
        print(f"  Type Capacity ({precision}): ~1e{limit}")
        print(f"  Conclusion: Coefficients are too large to maintain signal precision.")
        print(f"  Fix: Use multistage decimation (e.g. M=4 then M={M//4}).")
        return False

    return True


def design_polyphase_butterworth_numpy(order, cutoff, M):
    """
    Reference implementation.
    """
    z, p, k = signal.butter(order, cutoff, btype="low", output="zpk", analog=False)

    p_new = p**M
    denom_low = np.real(np.poly(p_new))

    s = np.array([1.0], dtype=complex)
    for pi in p:
        S_i = pi ** np.arange(M)
        s = np.convolve(s, S_i)
    s = np.real(s)  # Imaginary parts cancel out since poles occur in conjugate pairs

    b = np.real(np.poly(z)) * k
    q = np.convolve(b, s)

    polyphase_branches = []
    for k in range(M):
        coeffs = q[k::M]
        polyphase_branches.append(coeffs)

    return polyphase_branches, denom_low


def test_and_plot_decomposition(order=3, cutoff=0.5, M=3):
    """
    Tests the polyphase decomposition by reconstructing the transfer function
    and comparing it against the original.
    """
    print(f"--- Running Test (Order={order}, M={M}, Cutoff={cutoff}) ---")

    b_orig, a_orig = signal.butter(order, cutoff, btype="low")
    poly_branches, denom_low = design_polyphase_butterworth(order, cutoff, M)

    q_len = sum(len(br) for br in poly_branches)
    q_recon = np.zeros(q_len)
    for k, branch in enumerate(poly_branches):
        q_recon[k::M] = branch
    q_recon = np.trim_zeros(q_recon, "b")

    a_recon = np.zeros(len(denom_low) * M - (M - 1))
    a_recon[::M] = denom_low

    lhs = np.convolve(b_orig, a_recon)
    rhs = np.convolve(q_recon, a_orig)

    max_len = max(len(lhs), len(rhs))
    lhs_pad = np.pad(lhs, (0, max_len - len(lhs)))
    rhs_pad = np.pad(rhs, (0, max_len - len(rhs)))

    max_error = np.max(np.abs(lhs_pad - rhs_pad))
    print(f"Numerical Verification (Cross-Multiplication Error): {max_error:.2e}")

    if max_error < 1e-10:
        print(">> SUCCESS: Decomposition matches original filter.")
    else:
        print(">> FAILURE: Reconstruction mismatch.")

    w, h_orig = signal.freqz(b_orig, a_orig, worN=2048)
    w, h_recon = signal.freqz(q_recon, a_recon, worN=2048)
    freq_axis = w / np.pi

    plt.figure(figsize=(10, 8))

    def toDecibel(h_abs):
        return 20 * np.log10(h_abs)

    plt.subplot(2, 1, 1)
    plt.plot(
        freq_axis,
        toDecibel(np.abs(h_orig)),
        "b",
        linewidth=3,
        alpha=0.6,
        label="Original",
    )
    plt.plot(
        freq_axis,
        toDecibel(np.abs(h_recon)),
        "r--",
        linewidth=2,
        label="Reconstructed",
    )
    plt.axvline(cutoff, color="black", ls="--", alpha=0.5, label="cutoff")
    plt.title(f"Filter Response Comparison (Order {order}, M={M})")
    plt.ylabel("Magnitude (dB)")
    plt.ylim([-60, 5])

    plt.legend()
    plt.grid(True)

    plt.subplot(2, 1, 2)
    plt.plot(freq_axis, np.angle(h_orig), "b", linewidth=3, alpha=0.6, label="Original")
    plt.plot(freq_axis, np.angle(h_recon), "r--", linewidth=2, label="Reconstructed")
    plt.axvline(cutoff, color="black", ls="--", alpha=0.5, label="cutoff")
    plt.ylabel("Phase (radians)")
    plt.xlabel("Normalized Frequency (×π rad/T)")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()


def pad_branches_to_equal_length(branches, pad_value=0.0):
    """Pads a list of sequence branches to have uniform length."""
    max_len = max(len(b) for b in branches)
    padded = []
    for b in branches:
        padded_row = list(b) + [pad_value] * (max_len - len(b))
        padded.append(padded_row)
    return padded, max_len


def format_cpp_array_1d(arr, type_cast="T"):
    """Formats a 1D sequence of numerical/mpmath values as a C++ std::array."""
    return ", ".join(f"{type_cast}({float(x):.17e})" for x in arr)


def format_cpp_array_2d(matrix, type_cast="T", indent="    "):
    """Formats a 2D sequence of numerical/mpmath values as nested C++ arrays."""
    lines = []
    for row in matrix:
        formatted_row = ", ".join(f"{type_cast}({float(x):.17e})" for x in row)
        lines.append(f"{indent}{{ {formatted_row} }},")
    return "\n".join(lines)


def generate_polyphase_cpp_struct(
    branches, denom_sections, struct_name="HybridCoefficients"
):
    """Generates the C++ template struct code for PolyphaseIir coefficients."""
    padded_branches, max_len = pad_branches_to_equal_length(branches, 0.0)
    M = len(branches)
    denom_len = len(denom_sections)

    denom_str = format_cpp_array_2d(denom_sections, type_cast="T", indent="    ")
    branches_str = format_cpp_array_2d(padded_branches, type_cast="T", indent="    ")

    code = f"""template<typename T> struct {struct_name} {{
  static constexpr int nPhase = {M};

  // SOS Denominator: {{a1, a2}} for 1 / (1 + a1*z^-1 + a2*z^-2)
  static constexpr std::array<std::array<T, 2>, {denom_len}> denom{{{{
{denom_str}
  }}}};

  static constexpr std::array<std::array<T, {max_len}>, nPhase> branches{{{{
{branches_str}
  }}}};
}};"""
    return code


def generate_sos_cpp_struct(sos_matrix, struct_name="SosTest"):
    """Generates standard SOS cascade representation for reference."""
    sections = []
    for section in sos_matrix:
        b0, b1, b2, a0, a1, a2 = section
        # Normalize coefficients by a0
        coeffs = [b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0]
        sections.append(coeffs)

    sos_str = format_cpp_array_2d(sections, type_cast="T", indent="    ")
    code = f"""template<typename T> struct {struct_name} {{
  // SOS coefficients as [b0, b1, b2, a1, a2]
  static constexpr std::array<std::array<T, 5>, {len(sos_matrix)}> co{{{{
{sos_str}
  }}}};
}};"""
    return code


def generate_polyphase_sos_cpp_struct(sos_polyphase, struct_name="SosCoefficients"):
    M, NSOS = len(sos_polyphase), len(sos_polyphase[0])
    lines = []
    indent = "    "
    for p_idx, phase in enumerate(sos_polyphase):
        lines.append(indent + "{{ // Phase " + str(p_idx))
        for s_idx, sec in enumerate(phase):
            b0, b1, b2, a0, a1, a2 = sec
            b0_n, b1_n, b2_n, a1_n, a2_n = b0 / a0, b1 / a0, b2 / a0, a1 / a0, a2 / a0
            formatted_sec = ", ".join(
                "T(" + f"{float(x):.17e}" + ")" for x in [b0_n, b1_n, b2_n, a1_n, a2_n]
            )
            lines.append(indent + indent + "{{" + formatted_sec + "}},")
        lines.append(indent + "}},")

    code = "template<typename T> struct " + struct_name + " {\n"
    code += "  static constexpr int nPhase = " + str(M) + ";\n"
    code += "  static constexpr int nSections = " + str(NSOS) + ";\n\n"
    code += "  static constexpr std::array<std::array<std::array<T, 5>, nSections>, nPhase> co{{\n"
    code += "\n".join(lines) + "\n"
    code += "  }};\n"
    code += "};"
    return code


def generate_polyphase_sos2_cpp_struct(sos_polyphase, struct_name="SosCoefficients2"):
    M, NSOS = len(sos_polyphase), len(sos_polyphase[0])

    b_lines = []
    indent = "    "
    for p_idx, phase in enumerate(sos_polyphase):
        b_lines.append(indent + "{{ // Phase " + str(p_idx))
        for s_idx, sec in enumerate(phase):
            b0, b1, b2, a0, a1, a2 = sec
            b0_n, b1_n, b2_n = b0 / a0, b1 / a0, b2 / a0
            formatted_sec = ", ".join(
                "T(" + f"{float(x):.17e}" + ")" for x in [b0_n, b1_n, b2_n]
            )
            b_lines.append(indent + indent + "{{" + formatted_sec + "}},")
        b_lines.append(indent + "}},")

    a_lines = []
    for s_idx in range(NSOS):
        b0, b1, b2, a0, a1, a2 = sos_polyphase[0][s_idx]
        a1_n, a2_n = a1 / a0, a2 / a0
        formatted_sec = ", ".join("T(" + f"{float(x):.17e}" + ")" for x in [a1_n, a2_n])
        a_lines.append(indent + "{{" + formatted_sec + "}},")

    code = "template<typename T> struct " + struct_name + " {\n"
    code += "  static constexpr int nPhase = " + str(M) + ";\n"
    code += "  static constexpr int nSections = " + str(NSOS) + ";\n\n"
    code += "  // Numerators: b0, b1, b2 per phase and section\n"
    code += "  static constexpr std::array<std::array<std::array<T, 3>, nSections>, nPhase> b{{\n"
    code += "\n".join(b_lines) + "\n"
    code += "  }};\n\n"
    code += "  // Denominators: a1, a2 per section, shared across all phases\n"
    code += "  static constexpr std::array<std::array<T, 2>, nSections> a{{\n"
    code += "\n".join(a_lines) + "\n"
    code += "  }};\n"
    code += "};"
    return code


def print_polyphase_cpp_struct(order, cutoff, M, struct_name="ButterworthCoefficients"):
    branches, sos_sections = design_polyphase_butterworth(order, cutoff, M, "hybrid")
    print(generate_polyphase_cpp_struct(branches, sos_sections, struct_name))


def print_sos_cpp_struct_f64(order, cutoff, struct_name="SosTest"):
    sos_matrix = signal.butter(order, cutoff, btype="low", output="sos")
    print(generate_sos_cpp_struct(sos_matrix, struct_name))


def print_coefficient_counts_table(max_order=16):
    import sympy

    def count_elements(data):
        """Recursively count all non-list/non-tuple elements."""
        n = 0
        stack = [data]
        while stack:
            current = stack.pop()
            if isinstance(current, (list, tuple)):
                for item in current:
                    stack.append(item)
            else:
                n += 1
        return n

    m_range = list(range(1, max_order + 1))
    fc = 0.125

    hybrid_results = {}
    sos_results = {}
    sos2_results = {}

    for order in range(1, max_order + 1):
        for M in sympy.divisors(order):
            if M == 1:
                continue

            try:
                q, sos_a = design_polyphase_butterworth(
                    order, fc, M, output="hybrid", workdps=16
                )
            except ValueError:
                pass
            else:
                hybrid_results[(order, M)] = count_elements(q) + count_elements(sos_a)

            try:
                sos_full = design_polyphase_butterworth(
                    order, fc, M, output="sos", workdps=16
                )
            except ValueError:
                pass
            else:
                sos_results[(order, M)] = (count_elements(sos_full) // 6) * 5

            try:
                sos2b, sos2a = design_polyphase_butterworth(
                    order, fc, M, output="sos2", workdps=16
                )
            except ValueError:
                pass
            else:
                sos2_results[(order, M)] = count_elements(sos2b) + count_elements(sos2a)

    def render_markdown_table(title, results):
        lines = []
        lines.append(f"### {title}")
        lines.append("")
        # Header
        header = "Order \\ M | " + " | ".join(map(str, m_range))
        separator = "-:|" + "-:| " * len(m_range)
        separator = separator[:-2]
        lines.append(header)
        lines.append(separator)

        for order in range(1, max_order + 1):
            row = [str(order)]
            for M in m_range:
                val = results.get((order, M), "")
                row.append(str(val))
            lines.append(" | ".join(row))

        return "\n".join(lines)

    print(
        render_markdown_table(
            "Hybrid Coefficient Counts (q + sparse denominator)", hybrid_results
        )
    )
    print("\n")
    print(
        render_markdown_table(
            "SOS Coefficient Counts (Full SOS per branch)", sos_results
        )
    )
    print("\n")
    print(
        render_markdown_table(
            "SOS2 Coefficient Counts (Shared denominators)", sos2_results
        )
    )


if __name__ == "__main__":
    order = 16
    cutoff = 0.015625
    fold = 8
    # order = 8
    # cutoff = 0.125
    # fold = 4

    # test_and_plot_decomposition(order, cutoff, fold)
    # print_polyphase_cpp_struct(order, cutoff, fold)
    # print_sos_cpp_struct(order, cutoff)
    print_coefficient_counts_table(16)
