import numpy as np
import numpy.polynomial as polynomial
import scipy.optimize as optimize
import scipy.special as special
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
import sympy

from numpy.polynomial import Polynomial
from numpy.polynomial import Chebyshev


def poly3(x, a0, a1, a2, a3):
    return a0 + a1 * x + a2 * x * x + a3 * x * x * x


def rat2(x, b1, b2, a1, a2):
    return (b1 * x + b2 * x * x) / (1 + a1 * x + a2 * x * x)


def rat3(x, b1, b2, b3, a1, a2, a3):
    return (b1 * x + b2 * x * x + b3 * x * x * x) / (
        1 + a1 * x + a2 * x * x + a3 * x * x * x
    )


def rat4(x, b1, b2, b3, b4, a1, a2, a3, a4):
    return (b1 * x + b2 * x * x + b3 * x * x * x + b4 * x * x * x * x) / (
        1 + a1 * x + a2 * x * x + a3 * x * x * x + a4 * x * x * x * x
    )


def curvefit():
    xdata = np.linspace(0, 1, 1024)
    ydata = np.sin(0.5 * np.pi * xdata)

    sigma = np.full_like(xdata, 6)
    sigma[-1] = 1

    popt_rat2, _ = optimize.curve_fit(
        rat2,
        xdata,
        ydata,
        sigma=sigma,
        absolute_sigma=True,
    )
    popt_rat3, _ = optimize.curve_fit(
        rat3,
        xdata,
        ydata,
        sigma=sigma,
        absolute_sigma=True,
    )
    popt_rat4, _ = optimize.curve_fit(
        rat4,
        xdata,
        ydata,
        sigma=sigma,
        absolute_sigma=True,
    )

    plt.plot(xdata, ydata, label="target")
    # plt.plot(xdata, rat2(xdata, *popt_rat2), label="rat2")
    # plt.plot(xdata, rat3(xdata, *popt_rat3), label="rat3")
    plt.plot(xdata, rat4(xdata, *popt_rat4), label="rat4")
    plt.grid()
    plt.legend()
    plt.show()


def chebyfit():
    xdata = np.linspace(0, 1, 2**16 + 1)
    # ydata = np.sin((np.pi / 2) * xdata)
    ydata = np.cos((np.pi / 2) * xdata)

    coef = polynomial.chebyshev.chebfit(xdata, ydata, 7)
    cheby = Chebyshev(coef)
    approx_cheby = cheby(xdata)

    poly = cheby.convert(kind=Polynomial)
    coef = poly.coef
    coef[0] = 0
    coef /= Polynomial(coef)(1)
    if ydata[0] >= 1:
        coef = -coef  # Only used for cos.
        coef[0] = 1  # Only used for cos.
    poly = Polynomial(coef)
    approx_poly = poly(xdata)

    diff = ydata - approx_poly
    max_error = np.max(np.abs(diff))
    print(
        f"""chebyfit
  a   : {poly.coef}
  diff: {max_error}
  bits: {int(-np.log2(max_error))}
"""
    )
    # plt.plot(xdata, diff, label="diff")

    plt.plot(xdata, ydata, label="target")
    plt.plot(xdata, approx_cheby, label="chebyshev")
    plt.plot(xdata, approx_poly, label="polynomial")
    plt.grid()
    plt.legend()
    plt.show()


def chebyfitShifted():
    xdata = np.linspace(-1 / 2, 1 / 2, 2**16 + 1)
    ydata_cos = np.cos((np.pi / 2) * (xdata + 1 / 2))
    ydata_sin = np.sin((np.pi / 2) * (xdata + 1 / 2))

    def approx(name, ydata):
        coef = polynomial.chebyshev.chebfit(xdata, ydata, 7)
        cheby = Chebyshev(coef)
        # approx_cheby = cheby(xdata)

        poly = cheby.convert(kind=Polynomial)
        # coef = poly.coef
        # coef[0] =
        # coef /= Polynomial(coef)(1)
        # poly = Polynomial(coef)
        approx_poly = poly(xdata)

        print(poly(-1 / 2), poly(1 / 2))

        diff = ydata - approx_poly
        max_error = np.max(np.abs(diff))
        print(
            f"""{name}
  a   : {poly.coef}
  diff: {max_error}
  bits: {int(-np.log2(max_error))}
"""
        )
        return approx_poly

    approx_cos = approx("chebyfit shifted cos", ydata_cos)
    approx_sin = approx("chebyfit shifted sin", ydata_sin)

    # plt.plot(xdata, diff, label="diff")

    plt.plot(xdata, ydata_cos, label="target cos")
    plt.plot(xdata, ydata_sin, label="target sin")
    plt.plot(xdata, approx_cos, label="approx. cos")
    plt.plot(xdata, approx_sin, label="approx. sin")
    plt.grid()
    plt.legend()
    plt.show()


def pade():
    sin_taylor = [0]
    for k in range(1, 16):
        if k % 2 == 0:
            sin_taylor.append(0)
            continue
        sin_taylor.append((-1) ** (k - 1) / special.factorial(k))
    print(sin_taylor)

    p, q = interpolate.pade(sin_taylor, 5, 5)

    xdata = np.linspace(-2 * np.pi, 2 * np.pi, 1024)
    target = np.sin(xdata)
    value = p(xdata) / q(xdata)

    plt.plot(xdata, target, label="target")
    plt.plot(xdata, value, label="pade")
    plt.legend()
    plt.grid()
    plt.ylim((-2, 2))
    plt.show()


"""
Polynomial approximation.

## sin
y = a1 * x^1 + a3 * x^3 + a5 * x^5.

```
x2 = x * x;
y = x * (a1 + x2 * (a3 + x2 * a5));
```

A = [
      1,       1,       1, # point at 1.
  (1/2), (1/2)^3, (1/2)^5, # point at 1/2.
  (3/4), (3/4)^3, (3/4)^5, # point at 3/4.
]

b = [1, sin(pi/4), sin(pi*3/8)]^T

x = [a1, a3, a5]^T

Solve Ax = b for x. Something like `inv(A) * b`.

## cos
y = a0 + a2 * x^2 + a4 * x^4 + a6 * x^6.

```
x2 = x * x;
y = a0 + x2 * (a2 + x2 * (a4 + x2 * a6));
```
"""


def derivative(degree, n):
    with np.errstate(divide="ignore"):
        return np.where(
            degree - n < 0, 0, special.factorial(degree) / special.factorial(degree - n)
        )


def test_derivative():
    degrees = np.array([1, 3, 5])

    print(derivative(degrees, 1))
    print(derivative(degrees, 2))
    print(derivative(degrees, 3))


def sin_derivative_at_half_pi(n):
    return [1, 0, -1, 0][n % 4]


def sinApprox(points: np.array):
    """
    len(points) must be odd.
    """
    A = np.array([points**i for i in range(1, len(points) * 2, 2)]).T
    b = np.array([np.sin((np.pi / 2) * x) for x in points])
    return np.linalg.solve(A, b)


def sinApprox2(order=5):
    nPoints = int((order + 1) / 2)
    degrees = np.arange(1, order + 1, 2)

    A = [np.ones(nPoints)]
    b = [1]
    if order >= 2:
        A.append(derivative(degrees, 1))
        b.append(0)  # 1st derivative of sin at pi/2.

    index = 2
    while index < nPoints:
        x = (index - 1) / index
        A.append(x**degrees)
        b.append(np.sin((np.pi / 2) * x))
        index += 1

    return np.linalg.solve(A, b)


def plotQuadSin(coefficients):
    x = np.linspace(0, 1, 1024 + 1)
    target = np.sin((np.pi / 2) * x)

    ## `a` is coefficient of odd order polynomial.
    ## y = a[0] * x^1 + a[1] * x^3 + a[2] * x^5 + ...
    for data in coefficients:
        a = data["coeff"]
        name = data["name"]

        # derivatives at pi/2.
        deriv = []
        for d in range(1, 3):
            x0 = 1
            y0 = 0
            for i, a_n in enumerate(a, 1):
                degree = i * 2 - 1
                y0 += a_n * derivative(degree, d) * x0**degree
            deriv.append(y0)

        y = np.zeros_like(x)
        for i, a_n in enumerate(a, 1):
            y += a_n * x ** (i * 2 - 1)

        diff = target - y
        max_error = np.max(np.abs(diff))
        print(
            f"""{name}
  a       : {list(a)}
  d1 at 1 : {deriv}
  max diff: {max_error}
  acc. bit: {int(-np.log2(max_error))}
"""
        )
        plt.title("Error (target - approx.)")
        plt.plot(diff, label=name)
    #     plt.plot(y, label=name)
    # plt.plot(target, label="target")
    plt.grid()
    plt.legend()
    plt.show()


def createCoeff(name, a):
    return {"name": name, "coeff": a}


def testSin():
    order = 13
    N = int((order + 1) / 2)

    coefficients = [
        createCoeff(
            "(n-1)/n",
            sinApprox(np.array([1] + [(i - 1) / i for i in range(2, N + 1)])),
        ),
        # createCoeff(
        #     "1-1/2^n",
        #     sinApprox(np.array([1] + [1 - 0.5**i for i in range(2, N + 1)])),
        # ),
        createCoeff(
            "n/N",
            sinApprox(np.array(np.linspace(1 / N, 1, N))),
        ),
        createCoeff(
            "1st deriv.",
            sinApprox2(order),
        ),
    ]
    plotQuadSin(coefficients)


def cosApprox(points: np.array):
    """
    len(points) must be even.
    """
    A = np.array([points**i for i in range(0, len(points) * 2, 2)]).T
    b = np.array([np.cos((np.pi / 2) * x) for x in points])
    coeff = np.linalg.solve(A, b)
    return coeff / coeff[0]


def plotQuadCos(coefficients):
    x = np.linspace(0, 1, 1024 + 1)
    target = np.cos((np.pi / 2) * x)

    ## `a` is coefficient of odd order polynomial.
    ## y = a[0] * x^1 + a[1] * x^3 + a[2] * x^5 + ...
    for data in coefficients:
        a = data["coeff"]
        name = data["name"]

        # derivatives at pi/2.
        deriv = []
        for d in range(0, 3):
            x0 = 1
            y0 = 0
            for i, a_n in enumerate(a):
                degree = i * 2
                y0 += a_n * derivative(degree, d) * x0**degree
            deriv.append(y0)

        y = np.zeros_like(x)
        for i, a_n in enumerate(a):
            y += a_n * x ** (i * 2)

        diff = target - y
        max_error = np.max(np.abs(diff))
        print(
            f"""{name}
  a       : {list(a)}
  dx at 1 : {deriv}
  max diff: {max_error}
  acc. bit: {int(-np.log2(max_error))}
"""
        )
        plt.plot(diff, label=name)
    #     plt.plot(y, label=name)
    # plt.plot(target, label="target")
    plt.grid()
    plt.legend()
    plt.show()


def testCos():
    order = 16
    N = int(order / 2) + 1

    coefficients = [
        # createCoeff(
        #     "(n-1)/n",
        #     cosApprox(np.array([1] + [(i - 1) / i for i in range(2, N + 1)])),
        # ),
        createCoeff(
            "n/N",
            cosApprox(np.array(np.linspace(1 / N, 1, N))),
        ),
    ]
    plotQuadCos(coefficients)


def sinCosApprox(order):
    """
    `order` should be even.
    """
    nPoints = order + 1
    poly = [Polynomial([0 for _ in range(i)] + [1]) for i in range(nPoints)]

    # np.polynomial.set_default_printstyle("ascii")
    # for p in poly:
    #     print(p)

    A = [[1] + [0 for _ in range(1, nPoints)]]
    b = [np.sin(np.pi / 4)]
    if order >= 2:
        A.append([p(-1 / 2) for p in poly])
        b.append(1)  # sin at pi/2.

    if order >= 3:
        A.append([p(1 / 2) for p in poly])
        b.append(0)  # sin at 0.

    x = sympy.symbols("x")
    index = 3
    while index < nPoints:
        n_th = int((index - 1) / 2)

        value = 1 / 2 if index % 2 == 0 else -1 / 2
        A.append([p.deriv(n_th)(value) for p in poly])

        diffed = sympy.diff(sympy.cos(x * (np.pi / 2)), x, n_th)
        evaled = float(diffed.subs(x, value + 1 / 2))
        if np.abs(evaled) < 1:
            evaled = 0
        b.append(evaled)

        index += 1

    np.set_printoptions(linewidth=120, suppress=True)
    return np.linalg.solve(A, b)


def testSinCos():
    order = 16

    x = np.linspace(-1 / 2, 1 / 2, 2**16 + 1)
    target_cos = np.cos((np.pi / 2) * x + np.pi / 4)
    target_sin = np.sin((np.pi / 2) * x + np.pi / 4)

    a = sinCosApprox(order)
    poly_cos = Polynomial(a)(x)
    diff_cos = target_cos - poly_cos
    max_diff_cos = np.max(np.abs(diff_cos))

    # Flip odd term sign.
    poly_sin = Polynomial(a * np.resize([1, -1], len(a)))(x)
    diff_sin = target_sin - poly_sin
    max_diff_sin = np.max(np.abs(diff_sin))

    print(
        f"""Olli sincos order {order}
  a           : {list(a)}
  a 31bit     : {list(np.floor(a*2**31))}
  max diff cos: {max_diff_cos}
  max diff sin: {max_diff_sin}
  acc. bit cos: {int(-np.log2(max_diff_cos))}
  acc. bit sin: {int(-np.log2(max_diff_sin))}
"""
    )

    # plt.plot(poly_cos, label="approx. cos", color="red")
    # plt.plot(poly_sin, label="approx. sin", color="red")
    # plt.plot(target_sin, label="numpy.sin", color="blue", alpha=0.5)
    # plt.plot(target_cos, label="numpy.cos", color="orange", alpha=0.5)

    plt.plot(diff_cos, label="diff. cos")
    plt.plot(diff_sin, label="diff. sin")

    plt.grid()
    plt.legend()
    plt.show()


def shiftedSinApprox(points: np.array, shift: float):
    px = points + shift
    A = np.array([px**i for i in range(0, len(points))]).T
    b = np.array(np.sin(2 * np.pi * points))
    return np.linalg.solve(A, b)


def shiftedCosApprox(points: np.array, shift: float):
    px = points + shift
    A = np.array([px**i for i in range(0, len(points))]).T
    b = np.array(np.cos(2 * np.pi * points))
    return np.linalg.solve(A, b)


def shiftedSinApprox2(d0func, d1func, d2func, order=9):
    """Use odd order."""
    nPoints = order // 2 + order % 2
    points_top = np.linspace(0, 1, nPoints, endpoint=True)
    power = np.arange(order + 1)
    A_top = np.array([pt**power for pt in points_top])

    d1 = np.array([derivative(degree, 1) for degree in range(order + 1)])
    points_bottom = np.linspace(0, 1, nPoints, endpoint=True)
    A_bottom = np.array([d1 * pt ** (power - 1) for pt in points_bottom])
    A_bottom = np.where(np.isfinite(A_bottom), A_bottom, 0)
    A = np.array(np.vstack((A_top, A_bottom)))
    b = np.array(
        np.hstack(
            (
                [d0func(pt) for pt in points_top],
                [d1func(pt) for pt in points_bottom],
            )
        )
    )

    # Add 2nd derivative of sin(2*pi*x) for x when order is even.
    # This didn't work well.
    if order % 2 == 0:
        point = 3 / 4
        d2 = np.array([derivative(degree, 2) for degree in range(order + 1)])
        row = d2 * point ** (power - 2)
        row = np.where(np.isfinite(row), row, 0)
        A = np.vstack((A, row))
        b = np.vstack((b, d2func(point)))

    return np.linalg.solve(A, b)


def testShiftedSin():
    def d0sin(x):
        return np.sin(2 * np.pi * x)

    def d1sin(x):
        return 2 * np.pi * np.cos(2 * np.pi * x)

    def d2sin(x):
        return 4 * np.pi * (np.cos(2 * np.pi * x) - np.pi * x * np.sin(2 * np.pi * x))

    def d0cos(x):
        return np.cos(2 * np.pi * x)

    def d1cos(x):
        return -2 * np.pi * np.sin(2 * np.pi * x)

    def d2cos(x):
        return -4 * np.pi * (np.sin(2 * np.pi * x) + np.pi * x * np.cos(2 * np.pi * x))

    def normalizeCoefficients(a):
        poly = Polynomial(a)
        y = [poly(x) for x in np.linspace(0, 1, 5, endpoint=True)]
        return a / np.max(np.abs(y))

    order = 13
    shift = 0
    px = np.linspace(0, 1, order, endpoint=True)

    x = np.linspace(0, 1, 2**16 + 1)
    target_sin = np.sin(2 * np.pi * x)
    target_cos = np.cos(2 * np.pi * x)

    a_sin = shiftedSinApprox(px, shift)
    # a_sin = shiftedSinApprox2(d0sin, d1sin, d2sin, order)
    a_sin = normalizeCoefficients(a_sin)
    psin = Polynomial(a_sin)
    poly_sin = psin(x + shift)
    diff_sin = target_sin - poly_sin
    max_diff_sin = np.max(np.abs(diff_sin))

    a_cos = shiftedCosApprox(px, shift)
    # a_cos = shiftedSinApprox2(d0cos, d1cos, d2cos, order)
    a_cos = normalizeCoefficients(a_cos)
    pcos = Polynomial(a_cos)
    poly_cos = pcos(x + shift)
    diff_cos = target_cos - poly_cos
    max_diff_cos = np.max(np.abs(diff_cos))

    print(
        f"""shifted sin order {len(a_sin)-1}
  a_sin           : {list(a_sin)}
  max diff sin: {max_diff_sin}
  acc. bit sin: {int(-np.log2(max_diff_sin))}
  sin at 0    : {psin(0)}
  sin at 1/2*π: {psin(1/4)}
  sin at π    : {psin(1/2)}
  sin at 3/2*π: {psin(3/4)}
  sin at 2*π  : {psin(1)}

shifted cos order {len(a_cos)-1}
  a_cos           : {list(a_cos)}
  max diff cos: {max_diff_cos}
  acc. bit cos: {int(-np.log2(max_diff_cos))}
  cos at 0    : {pcos(0)}
  cos at 1/2*π: {pcos(1/4)}
  cos at π    : {pcos(1/2)}
  cos at 3/2*π: {pcos(3/4)}
  cos at 2*π  : {pcos(1)}
"""
    )

    # plt.plot(x, poly_sin, label="approx. sin", color="red")
    # plt.plot(x, poly_cos, label="approx. cos", color="red")
    # plt.plot(x, target_sin, label="numpy.sin", color="blue", alpha=0.5)
    # plt.plot(x, target_cos, label="numpy.cos", color="orange", alpha=0.5)

    plt.plot(x, diff_sin, label="diff. sin")
    plt.plot(x, diff_cos, label="diff. cos")

    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # curvefit()
    # chebyfit()
    # chebyfitShifted()
    # testSin()
    # testCos()
    # testSinCos()
    testShiftedSin()
