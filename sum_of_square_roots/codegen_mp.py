import math
import mpmath


def generate_sum_sqrt_code(dps: int) -> str:
    """Generates Python code for approximating the sum of square roots up to N

    with the specified decimal precision (dps).

    Parameters:
    -----------
    dps : int
        The target decimal precision.

    Returns:
    --------
    str
        The Python code defining the lazy-cached approximation function.
    """
    mpmath.mp.dps = dps + 10
    epsilon = mpmath.power(10, -dps)

    max_k = max(200, dps + 50)
    coeffs = []
    for k in range(1, max_k + 1):
        if k == 1:
            deriv = mpmath.mpf(0.5)
        else:
            deriv = coeffs[-1]["deriv"] * (0.5 - (2 * k - 3)) * (0.5 - (2 * k - 2))

        b_val = mpmath.bernoulli(2 * k)
        fact = mpmath.factorial(2 * k)
        total_coeff = b_val / fact * deriv
        coeffs.append({"coeff": total_coeff, "deriv": deriv})

    max_M_limit = 1000
    best_k = None
    best_M = None

    for k in range(1, max_k):
        power = 2 * k - 1.5
        coeff_val = coeffs[k - 1]["coeff"]

        M_prec = float((abs(coeff_val) / epsilon) ** (1.0 / power))

        if k > 1:
            M_conv = float(
                mpmath.sqrt(abs(coeffs[k - 1]["coeff"] / coeffs[k - 2]["coeff"]))
            )
        else:
            M_conv = 0.0

        M_req = max(M_prec, M_conv)
        if M_req <= max_M_limit:
            best_k = k
            best_M = max(1, int(math.ceil(M_req)))
            break

    if best_k is None:
        valid_pairs = []
        for k in range(1, max_k):
            power = 2 * k - 1.5
            coeff_val = coeffs[k - 1]["coeff"]
            M_prec = float((abs(coeff_val) / epsilon) ** (1.0 / power))
            M_conv = (
                float(mpmath.sqrt(abs(coeffs[k - 1]["coeff"] / coeffs[k - 2]["coeff"])))
                if k > 1
                else 0.0
            )
            M_req = max(M_prec, M_conv)
            valid_pairs.append((k, M_req))
        best_k, min_M_req = min(valid_pairs, key=lambda x: x[1])
        best_M = max(1, int(math.ceil(min_M_req)))

    selected_coeffs = [coeffs[i]["coeff"] for i in range(best_k)]
    zeta_val = mpmath.zeta(-0.5)

    zeta_str = mpmath.nstr(zeta_val, dps + 5)
    coeffs_str = [mpmath.nstr(c, dps + 5) for c in selected_coeffs]

    horner_lines = []
    if best_k > 0:
        horner_lines.append(f"    poly = mpmath.mpf('{coeffs_str[-1]}')")
        for c in reversed(coeffs_str[:-1]):
            horner_lines.append(f"    poly = poly * m + mpmath.mpf('{c}')")
    else:
        horner_lines.append("    poly = mpmath.mpf('0')")

    horner_code = "\n".join(horner_lines)

    code = f"""import mpmath

mpmath.mp.dps = {dps}

_sum_sqrt_cache = [mpmath.mpf(0)]

def sum_sqrt_approx(N):
    \"\"\"
    Approximation of sum(sqrt(i) for i in range(1, N+1)) using Euler-Maclaurin formula.
    Precision: {dps} dps
    Threshold M: {best_M}
    Order k: {best_k}
    \"\"\"
    if N < {best_M}:
        curr_len = len(_sum_sqrt_cache)
        if N >= curr_len:
            last_val = _sum_sqrt_cache[-1]
            for i in range(curr_len, N + 1):
                last_val = last_val + mpmath.sqrt(mpmath.mpf(i))
                _sum_sqrt_cache.append(last_val)
        return _sum_sqrt_cache[N]

    N_mp = mpmath.mpf(N)
    sqrt_n = mpmath.sqrt(N_mp)
    zeta_const = mpmath.mpf('{zeta_str}')

    base = (mpmath.mpf('2') / mpmath.mpf('3') * N_mp + mpmath.mpf('0.5')) * sqrt_n + zeta_const

    m = mpmath.mpf('1') / (N_mp * N_mp)

{horner_code}

    return base + (poly / sqrt_n)
"""
    return code


if __name__ == "__main__":
    print(generate_sum_sqrt_code(50))
