"""
zeta(-0.5) ~= -0.20788622497735457


f64: order 3 correction, LUT below N < 64.
f32: order 1 correction, LUT below N < 32.
"""

import math
from mpmath import mp
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
import sympy as sp

mp.dps = 720


def sum_sqrt_mp(N):
    return mp.nsum(lambda n: mp.sqrt(n), (1, N))


def sum_sqrt_mp_array(N):
    y = [mp.mpf(0)]
    for i in range(1, N):
        y.append(y[-1] + mp.sqrt(i))
    return y


def sum_sqrt_approx_mp(N, terms=1):
    zeta_term = mp.mpf(special.zeta(-0.5))

    N_mp = mp.mpf(N)
    res = (2 / 3) * mp.power(N_mp, 1.5) + 0.5 * mp.power(N_mp, 0.5) + zeta_term

    if terms >= 1:
        res += (1 / 24) * mp.power(N_mp, -0.5)
    if terms >= 2:
        res -= (1 / 1920) * mp.power(N_mp, -2.5)
    return res


def print_table(N, precision="f64"):
    if precision == "f32":
        convert = lambda x: np.float32(float(x))
        fmt = "1.8e"
    else:
        convert = float
        fmt = "1.17e"
    print("[" + ", ".join([f"{convert(sum_sqrt_mp(n)):{fmt}}" for n in range(N)]) + "]")


def print_sqrt_sum_coefficients(max_order=20):
    x = sp.Symbol("x")
    N = sp.Symbol("N")
    f = sp.sqrt(x)

    coeffs = []
    for k in range(2, max_order + 1, 2):
        P = sp.diff(f, x, k - 1).subs(x, 1)
        C = (sp.bernoulli(k) / sp.factorial(k)) * P
        power_of_N = sp.Rational(3, 2) - k
        coeffs.append({"order": k, "coefficient": C, "power_N": power_of_N})

    print(f"zeta(-0.5) ~= {mp.zeta(-0.5):+.17e}\n")
    print(
        f"{'Order (k)':<10} {'Power of N':<12} {'Coefficient (Rational)':<25} {'Coefficient (Float)'}"
    )
    print("-" * 75)

    expr = 0
    for item in coeffs:
        rat = item["coefficient"]
        rat_s = str(rat) if rat < 0 else " " + str(rat)  # Align the sign.
        pow_N = item["power_N"]

        print(f"{item['order']:<10} {str(pow_N): <12} {rat_s:<25} {float(rat):+.17e}")

        expr += rat * N**pow_N

    print(f"\n{sp.latex(expr)}")


def sum_sqrt_approx_order10(N):
    """
    Computes sum(sqrt(i) for i=1 to N) using Euler-Maclaurin approximation.
    Includes correction terms up to Order 10 (Bernoulli B22).

    Structure:
      S = Base_Terms + (1/sqrt(N)) * Polynomial(1/N^2)
    """
    if N < 6:
        return [
            0.0,
            1.0,
            2.414213562373095,
            4.146264369941973,
            6.146264369941973,
            8.382332347441762,
        ][N]

    sqrt_n = math.sqrt(N)
    zeta_const = -0.2078862249773545660  # zeta(-0.5)

    base = ((2.0 / 3.0) * N + (0.5)) * sqrt_n + zeta_const

    m = 1.0 / (N * N)

    poly = (-25273021529) / (274877906944)
    poly = poly * m + (4535189795) / (360777252864)
    poly = poly * m + (-4741887) / (2147483648)
    poly = poly * m + (52003) / (100663296)
    poly = poly * m + (-223193) / (1321205760)
    poly = poly * m + (65) / (786432)
    poly = poly * m + (-11) / (163840)
    poly = poly * m + (1) / (9216)
    poly = poly * m + (-1) / (1920)
    poly = poly * m + (1) / (24)

    return base + (poly / sqrt_n)


table64 = [
    0.00000000000000000e00,
    1.00000000000000000e00,
    2.41421356237309492e00,
    4.14626436994197256e00,
    6.14626436994197256e00,
    8.38233234744176237e00,
    1.08318220902249394e01,
    1.34775734012895310e01,
    1.63060005260357208e01,
    1.93060005260357208e01,
    2.24682781862040990e01,
    2.57849029765595006e01,
    2.92490045916972541e01,
    3.28545558671612454e01,
    3.65962132539351828e01,
    4.04691966001426024e01,
    4.44691966001426024e01,
    4.85923022257602639e01,
    5.28349429128795478e01,
    5.71938418564202209e01,
    6.16659778114198005e01,
    6.62485535063756430e01,
    7.09389692661990665e01,
    7.57348007895117945e01,
    8.06337802750781520e01,
    8.56337802750781520e01,
    9.07327997886709312e01,
    9.59289522113775632e01,
    1.01220454833506750e02,
    1.06605619640641251e02,
    1.12082845215692913e02,
    1.17650609578522932e02,
    1.23307463828015315e02,
    1.29052026474553344e02,
    1.34882978369398643e02,
    1.40799058152498247e02,
    1.46799058152498247e02,
    1.52881820682796473e02,
    1.59046234685765455e02,
    1.65291232684163845e02,
    1.71615788004500615e02,
    1.78018912241933464e02,
    1.84499652940341321e02,
    1.91057091464643321e02,
    1.97690341045354131e02,
    2.04398544977853476e02,
    2.11180874960978770e02,
    2.18036529561379808e02,
    2.24964732791655308e02,
    2.31964732791655308e02,
    2.39035800603520784e02,
    2.46177229032063622e02,
    2.53388331582991611e02,
    2.60668441472272150e02,
    2.68016910700621679e02,
    2.75433109187717321e02,
    2.82916423961265195e02,
    2.90466258396535977e02,
    2.98082031502399843e02,
    3.05763177250268484e02,
    3.13509143942683295e02,
    3.21319393618589970e02,
    3.29193401492601765e02,
    3.37130655425795567e02,
]


def sum_sqrt_approx_f64(N, dtype=np.float64):
    if N < 64:
        return dtype(table64[N])

    sqrt_n = np.sqrt(N, dtype=dtype)
    zeta_const = dtype(-0.2078862249773545660)  # zeta(-0.5)

    base = (dtype(2.0 / 3.0) * N + dtype(0.5)) * sqrt_n + zeta_const

    m = dtype(1) / (N * N)

    poly = dtype(1) / dtype(9216)
    poly = poly * m + dtype(-1) / dtype(1920)
    poly = poly * m + dtype(1) / dtype(24)

    return base + (poly / sqrt_n)


def sum_sqrt_approx_f32(N, dtype=np.float32):
    if N < 64:
        return dtype(table64[N])

    sqrt_n = np.sqrt(N, dtype=dtype)
    zeta_const = dtype(-0.2078862249773545660)  # zeta(-0.5)

    base = (dtype(2) / dtype(3) * N + dtype(0.5)) * sqrt_n + zeta_const
    return base + (dtype(1) / dtype(24) / sqrt_n)


def relative_error(ref, approx, dtype=np.float64):
    ref_arr = np.asarray(ref, dtype=dtype)
    approx_arr = np.asarray(approx, dtype=dtype)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = np.abs((ref_arr - approx_arr) / ref_arr, dtype=dtype)
    return np.nan_to_num(error, nan=0.0)


def ulp_error(ref, approx, dtype=np.float64):
    approx = np.asanyarray(approx, dtype=dtype)
    ref = np.asanyarray(ref, dtype=dtype)

    eps = np.spacing(ref, dtype=dtype)
    with np.errstate(divide="ignore", invalid="ignore"):
        error = np.abs(approx - ref, dtype=dtype) / eps

    error = np.where(approx == ref, 0.0, error)
    error = np.where(np.isinf(ref) & (approx != ref), np.inf, error)

    return error.item() if error.ndim == 0 else error


def test_sum_sqrt_approx():
    N = 65536
    fp_type = np.float64
    typename = np.dtype(fp_type).name
    fn_approx = sum_sqrt_approx_f64 if fp_type == np.float64 else sum_sqrt_approx_f32

    eps = np.finfo(fp_type).eps

    # approx = np.array([sum_sqrt_approx_order10(n) for n in range(N)])
    approx = np.array([fn_approx(n, fp_type) for n in range(N)])
    exact = np.array([fp_type(v) for v in sum_sqrt_mp_array(N)])

    abs_e = np.abs(approx - exact, dtype=fp_type)
    rel_e = relative_error(exact, approx, dtype=fp_type)
    ulp_e = ulp_error(exact, approx, dtype=fp_type)

    x_N = np.arange(N)

    fig, ax = plt.subplots(2, 2, figsize=(14, 8))
    fig.suptitle(f"Sum of Square Roots Approximation Errors ({typename})")

    ax[0][0].set_title(r"$f(N) = \sum_{i=1}^N \sqrt{i}$")
    ax[0][0].plot(approx, label="Approx.", color="red", alpha=0.5)
    ax[0][0].plot(exact, label="Exact", color="blue", alpha=0.5)

    ax[1][1].set_title(f"Relative Error")
    ax[1][1].plot(x_N, rel_e, label=f"Rel. error", color="orange", lw=1)
    ax[1][1].axhline(
        eps, label=f"{fp_type.__name__} eps.", color="black", alpha=0.5, ls="--", lw=1
    )
    ax[1][1].set_yscale("log")

    ax[0][1].set_title(f"ULP Error")
    ax[0][1].plot(x_N, ulp_e, label=f"ULP error", color="orange", lw=1)

    ax[1][0].set_title(f"Absolute Error")
    ax[1][0].plot(x_N, abs_e, label=f"Abs. error", color="orange", lw=1)
    ax[1][0].plot(
        exact * eps, label="exact * eps", color="black", alpha=0.5, ls="--", lw=1
    )

    for row in ax:
        for axis in row:
            axis.legend()
            axis.grid(which="both", color="#f8f8f8")
            axis.set_xlabel("N")
    plt.tight_layout()
    plt.savefig(f"img/error_{typename}.svg")
    plt.show()

    # print(f"N = {N}")
    # print(f"Exact Sum:     {exact:.15f}")
    # print(f"Approximation: {approx:.15f}")
    # print(f"Absolute Error: {error:.5e}")
    # print(f"Relative Error: {error/exact:.5e}")


def test_large_N():
    N = 123456789
    eps = np.finfo(np.float64).eps

    approx = sum_sqrt_approx_f64(N)
    exact = float(sum_sqrt_mp(N))

    abs_error = abs(approx - exact)
    rel_error = relative_error(exact, approx)

    print(f"Exact     : {exact:1.17e}")
    print(f"Approx.   : {approx:1.17e}")
    print(f"Abs. Error: {abs_error:1.17e}")
    print(f"Rel. Error: {rel_error:1.17e}")
    print(f"Accurate? : {rel_error < eps}")


def compute_min_n_table(max_order: int = 20, dtype=np.float64):
    B = special.bernoulli(max_order + 2)

    P = [1.0]
    value = 1.0
    for i in range(1, max_order + 3):
        value *= 3.0 / 2.0 - i
        P.append(value)

    results = []
    for k in range(2, max_order + 2, 2):
        coeff = abs((B[k + 2] / special.factorial(k + 2)) * P[k + 1])
        power = k + 0.5
        min_n = (coeff / np.finfo(dtype).eps) ** (1.0 / power)
        results.append(
            {
                "k": k,
                "Term": f"b_{k}",
                "Power": power,
                "Coeff": coeff,
                "Min_N": min_n,
            }
        )

    def format_min_n(min_n):
        if min_n >= 1e7:
            return f"N > {min_n:1.1e}"
        return f"{int(np.ceil(min_n))}"

    def lookup_table_size(min_n, dtype):
        size = np.ceil(min_n) * dtype(1).nbytes

        si_prefix = ["", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi", "Ri", "Qi"]
        index = 0
        while size >= 1024 and index < len(si_prefix) - 1:
            size /= 1024
            index += 1

        return f"{int(size):>4} {si_prefix[index]}B"

    print(f"--- {np.dtype(dtype).name}")
    print(
        f"{'k':<4} {'Term':<6} {'Power':<6} {'Coeff':<14} {'Min N':<6} {'LUT Size':<10}"
    )
    for row in results:
        print(
            f"{row['k']:<4} "
            f"{row['Term']:<6} "
            f"{row['Power']:<6.1f} "
            f"{row['Coeff']:<14.7e} "
            f"{format_min_n(row['Min_N']):<6} "
            f"{lookup_table_size(row['Min_N'], dtype):<10}"
        )


if __name__ == "__main__":
    # compute_min_n_table()
    # print_sqrt_sum_coefficients()
    # print_table(64)
    test_sum_sqrt_approx()
    # test_large_N()
