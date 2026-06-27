import sympy as sp


def tf():
    z, f, Q, X = sp.symbols("z f Q X")
    HP, BP, LP = sp.symbols("HP BP LP")

    eq1 = sp.Eq(HP, X - z**-1 * LP - z**-1 * BP / Q)
    eq2 = sp.Eq(BP, z**-1 * BP + f * HP)
    eq3 = sp.Eq(LP, z**-1 * LP + f * BP)

    sol = sp.solve([eq1, eq2, eq3], (HP, BP, LP))

    H_HP = sp.simplify(sol[HP] / X)
    H_BP = sp.simplify(sol[BP] / X)
    H_LP = sp.simplify(sol[LP] / X)

    print("H_HP:")
    sp.pprint(H_HP)
    print("\nH_BP:")
    sp.pprint(H_BP)
    print("\nH_LP:")
    sp.pprint(H_LP)


def cutoff_upper_bound():
    f = sp.Symbol("f", real=True)
    eq = sp.Eq(4 * f, (2 - f**2) * (f + sp.sqrt(f**2 + 8)))
    print(sp.solve(eq, f))


def maxmally_flat_denominator():
    a1, a2, w, x = sp.symbols("a1 a2 w x", real=True)

    re = 1 + a1 * sp.cos(w) + a2 * sp.cos(2 * w)
    im = -(a1 * sp.sin(w) + a2 * sp.sin(2 * w))

    mag_sq = re**2 + im**2
    expanded_trig = sp.expand_trig(mag_sq)
    fully_expanded = sp.expand(expanded_trig)
    cos_only = sp.expand(fully_expanded.subs(sp.sin(w) ** 2, 1 - sp.cos(w) ** 2))
    px = sp.collect(cos_only.subs(sp.cos(w), x), x)

    print("Polynomial P(x) where x = cos(w):")
    sp.pprint(px)
    print("-" * 30)

    coeff_x2 = sp.Poly(px, x).coeff_monomial(x**2)
    coeff_x1 = sp.Poly(px, x).coeff_monomial(x)
    coeff_x0 = sp.Poly(px, x).coeff_monomial(1)

    print(f"Coefficient of x^2: {coeff_x2}")
    print(f"Coefficient of x^1: {coeff_x1}")
    print(f"Coefficient of x^0: {coeff_x0}")


def f_range():
    f = sp.Symbol("f", real=True)

    lhs = (sp.sqrt(f**2 + 8) + f) / 4
    rhs = 2 * f / (4 - f**2)

    domain = sp.Interval(0, 2, right_open=True)
    solution = sp.solvers.inequalities.solve_univariate_inequality(
        lhs >= rhs, f, relational=False, domain=domain
    )

    print(solution)


if __name__ == "__main__":
    # tf()
    # cutoff_upper_bound()
    # maxmally_flat_denominator()
    f_range()
