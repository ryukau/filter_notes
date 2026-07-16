import math
import sympy as sp


def div_diff(points, func):
    """Computes the divided difference of a function."""
    n = len(points) - 1
    if n == 0:
        return func(points[0])
    return sp.simplify(
        (div_diff(points[1:], func) - div_diff(points[:-1], func))
        / (points[-1] - points[0])
    )


def get_b_splines(t, x, n):
    if n == 0:
        return []
    b_splines = []
    for i in range(n):

        def func_eval(y_val):
            j = x.index(y_val)
            if j <= i:
                return 0
            else:
                return (y_val - t) ** (n - 1)

        val = n * div_diff(x[: n + 1], func_eval)
        b_splines.append(sp.simplify(val))
    return b_splines


def test1(N=3):
    """Note: Slow when N >= 4."""
    x = sp.symbols(f"x0:{N+1}", real=True)
    t = sp.symbols("t", real=True)
    F = sp.Function("F")

    for n in range(N + 1):
        lhs = math.factorial(n) * div_diff(list(x[: n + 1]), F)

        if n == 0:
            rhs = F(x[0])
        else:
            f = sp.diff(F(t), t, n)
            B_splines = get_b_splines(t, list(x), n)
            rhs = sum(
                sp.integrate(B_splines[i] * f, (t, x[i], x[i + 1])) for i in range(n)
            )

        print(f"test1 n={n} Equivalence:", sp.simplify(lhs - rhs) == 0)


def test2(N=3):
    """Note: Slow when N >= 4."""
    u = sp.Symbol("u")
    t = sp.Symbol("t")
    x = sp.symbols(f"x0:{N+1}", real=True)
    coefs = sp.symbols(f"c0:{N+1}")
    f_poly = sum(c * t**i for i, c in enumerate(coefs))

    F = [f_poly]
    for _ in range(1, N + 1):
        F.append(sp.integrate(F[-1], t))

    for n in range(N + 1):
        if n == 0:
            lhs = div_diff([x[0]], lambda val: F[0].subs(t, val))
            rhs = lhs
        else:
            lhs = math.factorial(n) * div_diff(
                list(x[: n + 1]), lambda val: F[n].subs(t, val)
            )
            B_splines = get_b_splines(t, list(x), n)

            rhs_terms = []
            for i in range(n):
                dt_step = x[i + 1] - x[i]
                t_sub = x[i] + u * dt_step
                B_u = B_splines[i].subs(t, t_sub)
                f_u = f_poly.subs(t, t_sub)
                rhs_terms.append(dt_step * sp.integrate(B_u * f_u, (u, 0, 1)))

            rhs = sum(rhs_terms)

        print(f"test2 n={n} Equivalence:", sp.simplify(lhs - rhs) == 0)


if __name__ == "__main__":
    test1(3)
    test2(3)
