import sympy as sp
import numpy as np


def analyze_taylor_and_thresholds():
    x = sp.Symbol("x", real=True, positive=True)  # x in [0, 0.5).
    omega = sp.Symbol("omega", real=True, positive=True)
    eps = sp.Symbol("eps", real=True, positive=True)

    sn = sp.sin(sp.pi * x)
    f_expr = 2 * sn / (sp.sqrt(sn**2 + 1) + sn)

    taylor_order = 8
    taylor_x = f_expr.series(x, 0, taylor_order)
    poly_x = taylor_x.removeO()
    poly_omega = poly_x.subs(x, omega / (2 * sp.pi)).expand()

    print("Taylor series of f(x) around 0:")
    print(taylor_x)
    print("\nTaylor series in terms of omega = 2*pi*x:")
    print(poly_omega)

    print("\n" + "=" * 50)
    print("Optimal Branching Thresholds (omega <= threshold)")
    print("Based on keeping relative truncation error <= machine epsilon (eps)")
    print("=" * 50)

    eps_64 = float(np.finfo(np.float64).eps)
    eps_32 = float(np.finfo(np.float32).eps)

    coefficients = [1, -1 / 2, 1 / 12, 1 / 24, -11 / 480, -1 / 720, 139 / 40320]
    for k in range(1, len(coefficients)):
        next_coeff = sp.Rational(coefficients[k])
        threshold_expr = (eps / sp.Abs(next_coeff)) ** (sp.Rational(1, k))

        val_64 = float(threshold_expr.subs(eps, eps_64))
        val_32 = float(threshold_expr.subs(eps, eps_32))

        print(f"n_term={k}:")
        print(f"  Formula: omega <= {threshold_expr}")
        print(f"  float64: omega < {val_64:.4e}")
        print(f"  float32: omega < {val_32:.4e}\n")


if __name__ == "__main__":
    analyze_taylor_and_thresholds()
