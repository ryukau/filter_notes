import sympy


def verifySoftclip2():
    x, h, a_1, a_2 = sympy.symbols("x, h, a_1, a_2", real=True)

    ###########################################################
    # K1_src = (
    #     a_1**2 / 2
    #     - a_1 * h
    #     + h * x
    #     + (
    #         a_1**3
    #         - 3 * a_1**2 * a_2
    #         + 3 * a_1 * a_2**2
    #         - 3 * a_2**2 * x
    #         + 3 * a_2 * x**2
    #         - x**3
    #     )
    #     / (-12 * a_1 + 12 * h)
    # )
    # K1_simplified = (
    #     +a_1 * (a_1 / 2 - h)
    #     + h * x
    #     + ((a_1 - a_2) ** 3 - (x - a_2) ** 3) / (12 * (h - a_1))
    # )
    # print(sympy.expand(K1_src) - sympy.expand(K1_simplified))

    ###########################################################
    # L1_src = (
    #     a_1**2 / 2
    #     - a_1 * h
    #     + h * x
    #     + (a_1**3 - 3 * a_1**2 * a_2 + 3 * a_1 * a_2**2 - a_2**3) / (-12 * a_1 + 12 * h)
    # )
    # L1_simplified = +a_1 * (a_1 / 2 - h) + h * x + (a_1 - a_2) ** 3 / (12 * (h - a_1))
    # print(sympy.expand(L1_src) - sympy.expand(L1_simplified))
    # print(sympy.latex(L1_simplified))

    ##########################################################
    # K2_src = (
    #     -(a_1**3) / 3
    #     + a_1**2 * h / 2
    #     + a_1**2 * x / 2
    #     - a_1 * h * x
    #     + h * x**2 / 2
    #     + (
    #         -3 * a_1**4
    #         + 8 * a_1**3 * a_2
    #         + 4 * a_1**3 * x
    #         - 6 * a_1**2 * a_2**2
    #         - 12 * a_1**2 * a_2 * x
    #         + 12 * a_1 * a_2**2 * x
    #         - 6 * a_2**2 * x**2
    #         + 4 * a_2 * x**3
    #         - x**4
    #     )
    #     / (-48 * a_1 + 48 * h)
    # )
    # K2_simplified = (a_1**2 * (3 * x - 2 * a_1) + 3 * (h * (x - a_1) ** 2)) / 6 + (
    #     -((-a_1 + x) ** 2)
    #     * ((x + a_1) ** 2 - 2 * (2 * a_2 * (a_1 - a_2 + x) - (a_1 - a_2) ** 2))
    #     / (48 * (h - a_1))
    # )
    # print(sympy.expand(K2_src) - sympy.expand(K2_simplified))
    # # print(sympy.latex(K2_simplified))

    ###########################################################
    # L2_src = (
    #     -(a_1**3) / 3
    #     + a_1**2 * h / 2
    #     + a_1**2 * x / 2
    #     - a_1 * h * x
    #     + h * x**2 / 2
    #     + (
    #         -3 * a_1**4
    #         + 8 * a_1**3 * a_2
    #         + 4 * a_1**3 * x
    #         - 6 * a_1**2 * a_2**2
    #         - 12 * a_1**2 * a_2 * x
    #         + 12 * a_1 * a_2**2 * x
    #         + a_2**4
    #         - 4 * a_2**3 * x
    #     )
    #     / (-48 * a_1 + 48 * h)
    # )
    # L2_simplified = (a_1**2 * (3 * x - 2 * a_1) + 3 * (h * (x - a_1) ** 2)) / 6 + (
    #     (a_1 - a_2) ** 3 * (4 * x - 3 * a_1 - a_2)
    # ) / (48 * (h - a_1))
    # print(sympy.expand(L2_src) - sympy.expand(L2_simplified))
    # print(sympy.latex(L2_simplified))


def verifySoftclipN():
    x, A, C, r_c, x_c, x_s, S = sympy.symbols("x, A, C,r_c, x_c, x_s, S", real=True)
    β = sympy.Symbol("β", real=True, positive=True)

    ###########################################################
    # K1_src = (
    #     A * (-r_c + x_c) ** (β + 1) / (β + 1)
    #     - A * (-x + x_c) ** (β + 1) / (β + 1)
    #     - C * r_c
    #     + C * x
    #     + r_c**2 / 2
    # )
    # K1_simplified = (
    #     (A / (β + 1)) * ((x_c - r_c) ** (β + 1) - (x_c - x) ** (β + 1))
    #     + C * (x - r_c)
    #     + r_c**2 / 2
    # )
    # # print(sympy.expand(K1_src) - sympy.expand(K1_simplified))
    # print(sympy.latex(K1_simplified, order="none"))

    ###########################################################
    # L1_src = (
    #     A * (-r_c + x_c) ** (β + 1) / (β + 1)
    #     - C * r_c
    #     + C * x_c
    #     + S * x**2 / 2
    #     - S * x_c**2 / 2
    #     + r_c**2 / 2
    #     + x * (A * (x_c - x_s) ** β + C - S * x_s)
    #     - x_c * (A * (x_c - x_s) ** β + C - S * x_s)
    # )
    # L1_simplified = (
    #     +(x - x_c) * (A * (x_c - x_s) ** β + C - S * x_s)
    #     + A * (x_c - r_c) ** (β + 1) / (β + 1)
    #     + C * (x_c - r_c)
    #     + (S / 2) * (x**2 - x_c**2)
    #     + r_c**2 / 2
    # )
    # # print(sympy.expand(L1_src) - sympy.expand(L1_simplified))
    # print(sympy.latex(L1_simplified, order="none"))

    ##########################################################
    # K2_src = (
    #     -A * r_c * (-r_c + x_c) ** (β + 1) / (β + 1)
    #     + A * x * (-r_c + x_c) ** (β + 1) / (β + 1)
    #     - A * (-r_c + x_c) ** (β + 2) / ((β + 1) * (β + 2))
    #     + A * (-x + x_c) ** (β + 2) / ((β + 1) * (β + 2))
    #     + C * r_c**2 / 2
    #     - C * r_c * x
    #     + C * x**2 / 2
    #     - r_c**3 / 3
    #     + r_c**2 * x / 2
    # )
    # K2_simplified = (
    #     A
    #     * (
    #         (x - r_c) * (x_c - r_c) ** (β + 1) / (β + 1)
    #         + ((-x + x_c) ** (β + 2) - (-r_c + x_c) ** (β + 2)) / ((β + 1) * (β + 2))
    #     )
    #     + (C / 2) * (x - r_c) ** 2
    #     + r_c**2 * (x / 2 - r_c / 3)
    # )
    # # print(sympy.expand(K2_src) - sympy.expand(K2_simplified))
    # print(sympy.latex(K2_simplified, order="none"))

    ###########################################################
    L2_src = (
        -A * r_c * (-r_c + x_c) ** (β + 1) / (β + 1)
        + A * x_c * (-r_c + x_c) ** (β + 1) / (β + 1)
        - A * (-r_c + x_c) ** (β + 2) / ((β + 1) * (β + 2))
        + C * r_c**2 / 2
        - C * r_c * x_c
        + C * x_c**2 / 2
        + S * x**3 / 6
        - S * x_c**3 / 6
        - r_c**3 / 3
        + r_c**2 * x_c / 2
        + x**2 * (A * (x_c - x_s) ** β / 2 + C / 2 - S * x_s / 2)
        + x
        * (
            -2 * A * r_c * (-r_c + x_c) ** β / (2 * β + 2)
            - 2 * A * x_c * β * (x_c - x_s) ** β / (2 * β + 2)
            + 2 * A * x_c * (-r_c + x_c) ** β / (2 * β + 2)
            - 2 * A * x_c * (x_c - x_s) ** β / (2 * β + 2)
            - 2 * C * r_c * β / (2 * β + 2)
            - 2 * C * r_c / (2 * β + 2)
            - S * x_c**2 * β / (2 * β + 2)
            - S * x_c**2 / (2 * β + 2)
            + 2 * S * x_c * x_s * β / (2 * β + 2)
            + 2 * S * x_c * x_s / (2 * β + 2)
            + r_c**2 * β / (2 * β + 2)
            + r_c**2 / (2 * β + 2)
        )
        - x_c**2 * (A * (x_c - x_s) ** β / 2 + C / 2 - S * x_s / 2)
        - x_c
        * (
            -2 * A * r_c * (-r_c + x_c) ** β / (2 * β + 2)
            - 2 * A * x_c * β * (x_c - x_s) ** β / (2 * β + 2)
            + 2 * A * x_c * (-r_c + x_c) ** β / (2 * β + 2)
            - 2 * A * x_c * (x_c - x_s) ** β / (2 * β + 2)
            - 2 * C * r_c * β / (2 * β + 2)
            - 2 * C * r_c / (2 * β + 2)
            - S * x_c**2 * β / (2 * β + 2)
            - S * x_c**2 / (2 * β + 2)
            + 2 * S * x_c * x_s * β / (2 * β + 2)
            + 2 * S * x_c * x_s / (2 * β + 2)
            + r_c**2 * β / (2 * β + 2)
            + r_c**2 / (2 * β + 2)
        )
    )
    L2_simplified = (
        +A * (x_c - r_c) ** (β + 2) / (β + 1) * (1 - 1 / (β + 2))
        + (C / 2) * (x_c - r_c) ** 2
        + (S / 6) * (x**3 - x_c**3)
        + r_c**2 * (x_c / 2 - r_c / 3)
        + (x - x_c)
        * (
            +A * ((x_c - r_c) ** (β + 1) / (β + 1) - x_c * (x_c - x_s) ** β)
            + S * x_c * (x_s - x_c / 2)
            - C * r_c
            + r_c**2 / 2
            + (x + x_c) / 2 * (A * (x_c - x_s) ** β + C - S * x_s)
        )
    )
    print(sympy.ratsimp(sympy.expand(L2_src) - sympy.expand(L2_simplified)))
    # print(sympy.latex(L2_simplified, order="none"))
    # print(sympy.factor(r_c**2 - 2 * r_c * x + x**2))


if __name__ == "__main__":
    # verifySoftclip2()
    verifySoftclipN()
