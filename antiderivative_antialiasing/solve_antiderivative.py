"""
Solve some memoryless nonlinearities for antiderivative antialiasing.

Some expression couldn't be solved by SymPy.
"""

import sympy
import functools


def printAntiderivativeSingle(name, J0, target):
    J1 = sympy.integrate(J0, target)
    J2 = sympy.integrate(J1, target)

    print(f"--- {name}")
    print(J0, end="\n\n")
    print(J1, end="\n\n")
    print(J2, end="\n\n")

    # latex = True
    # if latex:
    #     text = f"- {name}\n$$\n\\begin{{aligned}}\n"
    #     for order, expr in enumerate([J0, J1, J2]):
    #         expr = sympy.simplify(expr)
    #         text += f"J^{order} f_{{\\mathrm{{{name}}}}} &= {sympy.latex(expr)}\\\\\n"
    #     text += f"\\end{{aligned}}\n$$\n"
    #     print(text)


def printAntiderivative():
    x = sympy.Symbol("x", real=True)
    β = sympy.Symbol("β", real=True)
    pa = functools.partial(printAntiderivativeSingle, target=x)

    pa("tanh", sympy.tanh(x))
    pa("atan", sympy.atan(x))
    pa("algebraic", x / (1 + sympy.Abs(x)))
    pa("softplus", sympy.log(1 + sympy.exp(x)))
    pa("swish", x / (1 + sympy.exp(-β * x)))
    pa("sinrunge", sympy.sin(2 * sympy.pi * x) / (1 + 10 * x * x))
    pa("exppoly", x**β * sympy.exp(-sympy.Abs(x)))
    pa("cosdecay", (1 - sympy.cos(x)) / x)
    pa("test_li2", -sympy.log(1 - x) / x)


def printAntiderivativeHardClip():
    x = sympy.Symbol("x", real=True)

    J0_case0 = x
    J0_case1 = sympy.sign(x)

    J1_case0 = sympy.integrate(J0_case0, x)
    J1_case1 = sympy.integrate(J0_case1, x)
    J1_case1 += J1_case0.subs(x, 1) - J1_case1.subs(x, 1)  # 積分定数

    J2_case0 = sympy.integrate(J1_case0, x)
    J2_case1 = sympy.integrate(J1_case1, x)
    J2_case1 += J2_case0.subs(x, 1) - J2_case1.subs(x, 1)  # 積分定数

    for expr in [J0_case0, J0_case1, J1_case0, J1_case1, J2_case0, J2_case1]:
        sympy.pprint(expr)
        print()


def printAntiderivativeSoftClip2():
    x, h, a_1, a_2 = sympy.symbols("x, h, a_1, a_2", real=True)

    # J:       |x| < a1,
    # K: a1 <= |x| < a2,
    # L: a2 <= |x|.
    J0 = x
    K0 = h + (a_2 - x) ** 2 / (4 * (a_1 - h))
    L0 = h

    J1 = sympy.integrate(J0, x)
    K1 = sympy.integrate(K0, x)
    L1 = sympy.integrate(L0, x)

    K1 += J1.subs(x, a_1) - K1.subs(x, a_1)
    L1 += K1.subs(x, a_2) - L1.subs(x, a_2)

    J2 = sympy.integrate(J1, x)
    K2 = sympy.integrate(K1, x)
    L2 = sympy.integrate(L1, x)

    K2 += J2.subs(x, a_1) - K2.subs(x, a_1)
    L2 += K2.subs(x, a_2) - L2.subs(x, a_2)

    print("--- softclip2")
    for cases in [[J0, K0, L0], [J1, K1, L1], [J2, K2, L2]]:
        print("---")
        for expr in cases:
            expr = sympy.ratsimp(expr)
            # sympy.pprint(expr)
            # print(sympy.latex(expr))
            print(expr)


def printAntiderivativeSoftClipN():
    x, A, C, r_c, x_c, x_s, S = sympy.symbols("x, A, C,r_c, x_c, x_s, S", real=True)
    β = sympy.Symbol("β", real=True, positive=True)

    # J:        |x| <= r_c,
    # K: r_c <  |x| <  x_c,
    # L: x_c <= |x|.
    J0 = x
    K0 = C + A * (x_c - x) ** β
    L0 = C + A * (x_c - x_s) ** β + S * (x - x_s)

    J1 = sympy.integrate(J0, x)
    K1 = sympy.integrate(K0, x)
    L1 = sympy.integrate(L0, x)

    K1 += J1.subs(x, r_c) - K1.subs(x, r_c)
    L1 += K1.subs(x, x_c) - L1.subs(x, x_c)

    J2 = sympy.integrate(J1, x)
    K2 = sympy.integrate(K1, x)
    L2 = sympy.integrate(L1, x)

    K2 += J2.subs(x, r_c) - K2.subs(x, r_c)
    L2 += K2.subs(x, x_c) - L2.subs(x, x_c)

    print("--- softclipN")
    for cases in [[J0, K0, L0], [J1, K1, L1], [J2, K2, L2]]:
        print("---")
        for expr in cases:
            # expr = sympy.ratsimp(expr)
            # sympy.pprint(expr)
            # print(sympy.latex(expr))
            print(expr)


if __name__ == "__main__":
    printAntiderivative()
    printAntiderivativeHardClip()
    printAntiderivativeSoftClip2()
    printAntiderivativeSoftClipN()
