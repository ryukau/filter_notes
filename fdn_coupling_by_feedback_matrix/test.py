import numpy as np
from scipy.stats import ortho_group, special_ortho_group
from tqdm import tqdm
import sympy
import json


def testOrthogonal():
    """
    Test a random orthogonal matrix made from 2 smaller orthogonal matrices.
    """

    def check(matrix: np.ndarray, seed: int):
        try:
            np.testing.assert_almost_equal(
                matrix @ matrix.T, np.identity(matrix.shape[0])
            )
        except Exception as e:
            print(f"Seed: {seed}", e)

    for _ in tqdm(range(1000)):
        rngSeed = np.random.default_rng()

        seed = rngSeed.integers(0, np.iinfo(np.int64).max)
        rng = np.random.default_rng(seed)

        size = 16
        A = ortho_group.rvs(size, random_state=rng)
        D = ortho_group.rvs(size, random_state=rng)

        ratio = 1
        gain = 1 / np.sqrt(ratio * ratio + 1)

        F1 = gain * np.block([[ratio * A, D], [-A, ratio * D]])
        F2 = gain * np.block([[ratio * A, -D], [A, ratio * D]])

        check(F1, seed)
        check(F2, seed)


def testTranspose():
    """Test NumPy `@` operator."""
    v = np.arange(9)
    m1 = np.reshape(v, (3, 3))
    m2 = np.reshape(v * v, (3, 3))
    print(m1 @ m2.T)
    print(m2 @ m1.T)


def solveCoefficientsFull(dim: int, toString: bool = False) -> dict:
    """
    Solve for coefficients that can be used to construct an orthogonal matrix made from smaller orthogonal matrices.

    - `dim` is number of smaller matrices.
    - If `toString` is `True`, output is formatted for JSON output. Otherwise, the output is a dict of SymPy expressions.
    """
    print(dim)
    a = sympy.Symbol("a")
    b = sympy.Symbol("b")
    y = sympy.Symbol("y")

    # For general solution, change `h` to following.
    #
    # ```
    # h = {
    #     i: {j: sympy.Symbol(f"h_{i},{j}") if i != j else a for j in range(1, dim + 1)}
    #     for i in range(1, dim + 1)
    # }
    # ```
    #
    h = {
        i: {j: b if i != j else a for j in range(1, dim + 1)} for i in range(1, dim + 1)
    }

    # # If `dim == 3`, uncomment this block to get some answers.
    # h[1][dim] = b
    # h[dim][1] = b

    eq = []
    for i in range(1, dim + 1):
        for j in range(i, dim + 1):
            expr = 0
            for k in range(1, dim + 1):
                expr += h[i][k] * h[j][k]
            if i == j:
                expr -= y
            eq.append(expr)

    variables = set([elem for row in h.values() for elem in row.values()])
    variables.remove(a)
    variables.add(y)

    # result = sympy.solve(eq, *variables, dict=True)
    result = sympy.nonlinsolve(eq, list(variables))
    # print(sympy.latex(result)) # debug

    if toString:
        result = [{str(k): str(v) for k, v in r.items()} for r in result]

    return result


def solveCoefficientsMostlyDiag(dim: int, toString: bool = False) -> dict:
    """
    It seems that the solutions only exist when `dim` is 2 or odd (2, 3, 5, 7, ...).
    And the solution is all 0 except diagonal elements.
    """
    a = sympy.Symbol("a")
    y = sympy.Symbol("y")
    h = {
        i: {j: sympy.Symbol(f"h_{i},{j}") if i != j else a for j in range(1, dim + 1)}
        for i in range(1, dim + 1)
    }

    def confinedRange(dim, start, end):
        if start < 1:
            start = 1
        if end > dim + 1:
            end = dim + 1
        return range(start, end)

    eq = set([])
    for i in range(1, dim + 1):
        for j in range(i - 1, i + 2):
            if j < 1 or j > dim:
                continue

            expr = 0
            if j == i - 1:
                for k in confinedRange(dim, i - 1, i + 1):
                    expr += h[i][k] * h[j][k]
            elif j == i:
                for k in confinedRange(dim, i - 1, i + 2):
                    expr += h[i][k] * h[j][k]
            elif j == i + 1:
                for k in confinedRange(dim, i, i + 2):
                    expr += h[i][k] * h[j][k]
            else:
                raise Exception("Shouldn't reach here, but just in case.")

            if i == j:
                expr -= y

            eq.add(expr)

    variables = set([elem for row in h.values() for elem in row.values()])
    variables.remove(a)
    variables.add(y)

    result = sympy.solve(eq, *variables, dict=True)
    # print(sympy.latex(result))  # debug

    if toString:
        result = [{str(k): str(v) for k, v in r.items()} for r in result]

    return result


def testKronOrtho():
    def check(matrix: np.ndarray, seed: int):
        try:
            np.testing.assert_almost_equal(
                matrix @ matrix.T, np.identity(matrix.shape[0])
            )
        except Exception as e:
            print(f"Seed: {seed}", e)

    for _ in tqdm(range(1000)):
        rngSeed = np.random.default_rng()

        seed = 56729  #  rngSeed.integers(0, np.iinfo(np.int64).max)
        rng = np.random.default_rng(seed)

        size = 2
        A = ortho_group.rvs(size, random_state=rng)
        B = ortho_group.rvs(size, random_state=rng)
        C = np.kron(A, B)
        check(C, seed)

        print(C @ C.T)
        exit()


if __name__ == "__main__":
    # testOrthogonal()
    # testTranspose()
    print(solveCoefficientsFull(3, False))
    # print(solveCoefficientsMostlyDiag(3, False))
    # testKronOrtho()

    # # Write solutions to json.
    # results = {dim: solveCoefficientsFull(dim, True) for dim in range(2, 128 + 1)}
    # with open("out.json", "w", encoding="utf-8") as fp:
    #     json.dump(results, fp)
