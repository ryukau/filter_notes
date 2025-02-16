"""
Solve integrals used for modified version of antiderivative antialiasing.

In the following paper, triangular window was used. In this code, the window function was replaced to cosine window.

- Parker, J. D., Zavalishin, V., & Le Bivic, E. (2016, September). Reducing the aliasing of nonlinear waveshaping using continuous-time convolution. In Proc. Int. Conf. Digital Audio Effects (DAFx-16), Brno, Czech Republic (pp. 137-144).
"""

import sympy


def solveIllConditionCase():
    tau = sympy.symbols("τ", real=True)

    def M1(h):
        numer = sympy.integrate(tau * h, (tau, 0, 1))
        denom = sympy.integrate(h, (tau, 0, 1))
        print(sympy.latex(numer / denom))
        # print(sympy.latex(denom))

    M1((1 - sympy.cos(sympy.pi * tau)) / 2)
    M1((1 + sympy.cos(sympy.pi * tau)) / 2)

    # # Below is to verify the ill-conditioned case of `h_lin` in Parker's paper.
    # M1(tau)
    # M1(1 - tau)


def solveCosine_Tanh():
    tau, x0, x1 = sympy.symbols("τ, x_0, x_1")
    expr = (1 - sympy.cos(sympy.pi * tau)) / 2 * sympy.tanh(x0 + tau * (x1 - x0))
    result = sympy.integrate(expr, (tau, 0, 1))

    print(expr)
    print()
    print(result)


def solveCosine_HalfRect():
    tau, x0, x1, x2 = sympy.symbols("τ, x_0, x_1, x_2", real=True)

    print("--- Term 1")
    expr = (1 - sympy.cos(sympy.pi * tau)) / 2 * (x0 + tau * (x1 - x0))
    p = -x0 / (x1 - x0)
    cases = [
        0,
        sympy.integrate(expr, [tau, p, 1]),
        sympy.integrate(expr, [tau, 0, p]),
        sympy.integrate(expr, [tau, 0, 1]),
    ]
    for index, expr in enumerate(cases):
        expr = sympy.simplify(expr)
        # expr = sympy.latex(expr)
        print(f"-- Case {index}")
        print(expr, end="\n\n")

    print("--- Term 2")
    expr = (1 + sympy.cos(sympy.pi * tau)) / 2 * (x1 + tau * (x2 - x1))
    p = -x1 / (x2 - x1)
    cases = [
        0,
        sympy.integrate(expr, [tau, p, 1]),
        sympy.integrate(expr, [tau, 0, p]),
        sympy.integrate(expr, [tau, 0, 1]),
    ]
    for index, expr in enumerate(cases):
        expr = sympy.simplify(expr)
        # expr = sympy.latex(expr)
        print(f"-- Case {index}")
        print(expr, end="\n\n")


def solveCosine_Hardclip():
    tau, x0, x1, x2 = sympy.symbols("τ, x_0, x_1, x_2", real=True)

    print("--- Term 1")
    h_cos = (1 - sympy.cos(sympy.pi * tau)) / 2
    exprN = -h_cos
    exprP = h_cos
    exprL = h_cos * (x0 + tau * (x1 - x0))
    a = (-x0 - 1) / (x1 - x0)
    b = (-x0 + 1) / (x1 - x0)
    J = sympy.integrate
    cases = [
        -sympy.Rational(1, 2),
        J(exprN, [tau, 0, a]) + J(exprL, [tau, a, 1]),
        J(exprN, [tau, 0, a]) + J(exprL, [tau, a, b]) + J(exprP, [tau, b, 1]),
        J(exprL, [tau, 0, a]) + J(exprN, [tau, a, 1]),
        J(exprL, [tau, 0, 1]),
        J(exprL, [tau, 0, b]) + J(exprP, [tau, b, 1]),
        J(exprP, [tau, 0, b]) + J(exprL, [tau, b, a]) + J(exprN, [tau, a, 1]),
        J(exprP, [tau, 0, b]) + J(exprL, [tau, b, 1]),
        sympy.Rational(1, 2),
    ]
    for index, expr in enumerate(cases):
        print(f"-- Case {index}")
        expr = sympy.expand(expr)
        expr = sympy.cancel(expr)
        expr = sympy.trigsimp(expr)
        print(expr, end="\n\n")
        # print(sympy.latex(expr), end="\n\n")

    print("--- Term 2")
    h_cos = (1 + sympy.cos(sympy.pi * tau)) / 2
    exprN = -h_cos
    exprP = h_cos
    exprL = h_cos * (x1 + tau * (x2 - x1))
    a = (-x1 - 1) / (x2 - x1)
    b = (-x1 + 1) / (x2 - x1)
    J = sympy.integrate
    cases = [
        -sympy.Rational(1, 2),
        J(exprN, [tau, 0, a]) + J(exprL, [tau, a, 1]),
        J(exprN, [tau, 0, a]) + J(exprL, [tau, a, b]) + J(exprP, [tau, b, 1]),
        J(exprL, [tau, 0, a]) + J(exprN, [tau, a, 1]),
        J(exprL, [tau, 0, 1]),
        J(exprL, [tau, 0, b]) + J(exprP, [tau, b, 1]),
        J(exprP, [tau, 0, b]) + J(exprL, [tau, b, a]) + J(exprN, [tau, a, 1]),
        J(exprP, [tau, 0, b]) + J(exprL, [tau, b, 1]),
        sympy.Rational(1, 2),
    ]
    for index, expr in enumerate(cases):
        print(f"-- Case {index}")
        expr = sympy.expand(expr)
        expr = sympy.cancel(expr)
        expr = sympy.trigsimp(expr)
        print(expr, end="\n\n")
        # print(sympy.latex(expr), end="\n\n")


if __name__ == "__main__":
    # solveIllConditionCase()
    # solveCosine_Tanh()
    solveCosine_HalfRect()
    # solveCosine_Hardclip()
