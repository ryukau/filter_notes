import math
import sympy
from sympy.polys.orthopolys import chebyshevt_poly
from sympy.polys.polyfuncs import horner

cosineSumWindowCoefficients = {
    "blackman": [7938 / 18608, 9240 / 18608, 1430 / 18608],
    "nuttall": [0.355768, 0.487396, 0.144232, 0.012604],
    "blackmanharris": [0.35875, 0.48829, 0.14128, 0.01168],
    "blackmannuttall": [0.3635819, 0.4891775, 0.1365995, 0.0106411],
    "flattop": [0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947368],
}

u0 = sympy.symbols("u0")
n = sympy.symbols("n")
for key, a_ in cosineSumWindowCoefficients.items():
    expr = 0
    for i in range(len(a_)):
        expr += (-1) ** i * a_[i] * chebyshevt_poly(i, x=u0)
    print(f"{key} = {horner(expr)}")

    # Below is testing `expr` against the original form of cosine-sum function.
    N = 255
    direct = a_[0]
    for i in range(1, len(a_)):
        direct += (-1) ** i * a_[i] * sympy.cos(i * 2 * sympy.pi * n / N)

    for i in range(N + 1):
        poly = sympy.N(expr.subs({u0: sympy.cos(2 * sympy.pi * i / N)}))
        target = sympy.N(direct.subs({n: i}))
        if not math.isclose(poly, target, abs_tol=2**-53):
            print(f"Discrepancy at: n={i}, target={target}, poly={poly}")
