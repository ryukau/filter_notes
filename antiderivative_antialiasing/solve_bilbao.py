"""
This code prints out higher order expressions of antiderivative antialiasing (ADAA).

The expressions are mostly useless because real implemetation requires branching to avoid division by zero, and that's the problematic part of ADAA. Most of the time, it's better to use 1st order ADAA with a bit of oversampling, unless the nonlinearity introduces severe aliasing.

Reference: "Antiderivative Antialiasing for Memoryless Nonlinearities", Stefan Bilbao, Fabián Esqueda, Julian D. Parker, Vesa Välimäki, IEEE Signal Processing Letters, Vol. 24, No. 7, July 2017
"""

import sympy
import functools


def shiftP(expr, n):
    """`e_+` in reference, eq.10."""
    return expr.subs(n, n + 1)


def shiftN(expr, n):
    """`e_-` in reference, eq.10."""
    return expr.subs(n, n - 1)


def shiftM(expr, n, m):
    return expr.subs(n, n + m)


def diffP(expr, n):
    """`δ_+` in reference, eq.10."""
    return shiftP(expr, n) - expr


def diffN(expr, n):
    """`δ_-` in reference, eq.10."""
    return expr - shiftN(expr, n)


def diffOneSided(expr, n, x):
    """`D_-` in reference, eq.11."""
    return diffN(expr, n) / diffN(x[n], n)


def diffCentered(expr, n, x):
    """`e_- D_2` in reference, appears in eq.13."""
    numer = 2 * diffP(diffOneSided(expr, n, x), n)
    denom = shiftP(x[n], n) - shiftN(x[n], n)
    return numer / denom


def adaa(p, n, x, F):
    m = sympy.floor(p / 2)
    r = sympy.Mod(p, 2)
    expr = F[n]
    for _ in range(m):
        expr = diffCentered(expr, n, x)
    if r != 0:
        expr = diffOneSided(expr, n, x)
    expr = shiftM(expr, n, -m)
    return expr


def test():
    n = sympy.Idx("n")
    F = sympy.IndexedBase("F")
    x = sympy.IndexedBase("x")

    assert shiftP(F[n], n) == F[n + 1]
    assert shiftN(F[n], n) == F[n - 1]
    assert diffP(F[n], n) == F[n + 1] - F[n]
    assert diffN(F[n], n) == F[n] - F[n - 1]

    # Eq.14 in reference.
    assert diffOneSided(F[n], n, x) == (F[n] - F[n - 1]) / (x[n] - x[n - 1])

    # Eq.15 in reference.
    assert sympy.expand(shiftN(diffCentered(F[n], n, x), n)) == sympy.expand(
        (2 / (x[n] - x[n - 2]))
        * (
            (F[n] - F[n - 1]) / (x[n] - x[n - 1])
            - (F[n - 1] - F[n - 2]) / (x[n - 1] - x[n - 2])
        )
    )

    # Eq.16 in reference.
    assert sympy.expand(adaa(3, n, x, F)) == sympy.expand(
        (1 / (x[n - 1] - x[n - 2]))
        * (
            (2 / (x[n] - x[n - 2]))
            * (
                (F[n] - F[n - 1]) / (x[n] - x[n - 1])
                - (F[n - 1] - F[n - 2]) / (x[n - 1] - x[n - 2])
            )
            - (2 / (x[n - 1] - x[n - 3]))
            * (
                (F[n - 1] - F[n - 2]) / (x[n - 1] - x[n - 2])
                - (F[n - 2] - F[n - 3]) / (x[n - 2] - x[n - 3])
            )
        )
    )


if __name__ == "__main__":
    # test()

    # Print out lengthy equations.
    n = sympy.Idx("n")
    F = sympy.IndexedBase("F")
    x = sympy.IndexedBase("x")
    print(sympy.latex(adaa(1, n, x, F)), end="\n\n")
    print(sympy.latex(adaa(2, n, x, F)), end="\n\n")
    print(sympy.latex(adaa(3, n, x, F)), end="\n\n")
    print(sympy.latex(adaa(4, n, x, F)), end="\n\n")
    print(sympy.latex(adaa(5, n, x, F)), end="\n\n")
