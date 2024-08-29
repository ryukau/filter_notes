import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import numpy.polynomial as polynomial
from numpy.polynomial import Polynomial
from numpy.polynomial import Chebyshev
from dataclasses import dataclass


"""
```maxima
ω: %pi - acos(a / 2);
t: tan(ω / 2);
expr: (t - 1) / (t + 1);
```

```wolfram
(tan(π - acos(a / 2) / 2) - 1) / (tan(π - acos(a / 2) / 2) + 1)
```
"""


def targetFunc(a):
    """`a` is in [-2, 2]."""
    omega_a = np.pi - np.arccos(a / 2)
    t = np.tan(omega_a / 2)
    return (t - 1) / (t + 1)


def poly3(x, a1, a2, a3):
    return (a1 + (a2 + a3 * x * x) * x * x) * x


def poly4(x, a1, a2, a3, a4):
    xx = x * x
    return (a1 + (a2 + (a3 + a4 * xx) * xx) * xx) * x


def poly5(x, a1, a2, a3, a4, a5):
    xx = x * x
    return (a1 + (a2 + (a3 + (a4 + a5 * xx) * xx) * xx) * xx) * x


def rat2(x, b1, b2, a1, a2):
    return (b1 * x + b2 * x * x) / (1 + a1 * x + a2 * x * x)


def rat3(x, b1, b2, b3, a1, a2, a3):
    return (b1 * x + b2 * x * x + b3 * x * x * x) / (
        1 + a1 * x + a2 * x * x + a3 * x * x * x
    )


def rat4(x, b1, b2, b3, b4, a1, a2, a3, a4):
    return ((((b4 * x + b3) * x + b2) * x + b1) * x) / (
        (((a4 * x + a3) * x + a2) * x + a1) * x + 1
    )


def rat5(x, b1, b2, b3, b4, b5, a1, a2, a3, a4, a5):
    return (((((b5 * x + b4) * x + b3) * x + b2) * x + b1) * x) / (
        ((((a5 * x + a4) * x + a3) * x + a2) * x + a1) * x + 1
    )


@dataclass
class FitResult:
    label: str
    value: np.array


def impl(x):
    b = [0.24324073816096012, -0.2162512785673542, 0.0473854827531673]
    a = [-0.9569399371569998, 0.23899241566503945, -0.005191170339007643]
    return (x * (b[0] + x * (b[1] + x * b[2]))) / (
        1 + x * (a[0] + x * (a[1] + x * a[2]))
    )


def approx():
    xdata = np.linspace(0, 2, 1024)
    ydata = targetFunc(xdata)

    sigma = np.full_like(xdata, 6)
    sigma[-1] = 1

    def fit(fitFunc):
        popt, _ = optimize.curve_fit(fitFunc, xdata, ydata, sigma=sigma)
        print(list(popt))
        return FitResult(fitFunc.__name__, fitFunc(xdata, *popt))

    def fitCheby(order=9):
        coef = polynomial.chebyshev.chebfit(xdata, ydata, order)
        cheby = Chebyshev(coef)
        return cheby(xdata)

    data = [
        # FitResult("cheby", fitCheby()),
        # fit(rat2),
        fit(rat3),
        # fit(rat4), # Error is too large.
        # fit(rat5), # Can't coverge.
        # fit(poly3),
        # fit(poly4),
        # fit(poly5),
        FitResult("impl", impl(xdata)),
    ]

    fig, ax = plt.subplots(2)
    fig.set_size_inches((8, 8))
    ax[0].plot(xdata, ydata, color="black", label="target")
    ax[1].plot(xdata, ydata / ydata, color="black", label="target")

    cmap = plt.get_cmap("viridis")
    for index, dt in enumerate(data):
        label = dt.label
        value = dt.value
        color = cmap(index / len(data))
        ax[0].plot(xdata, value, color=color, alpha=0.75, label=label)
        ax[1].plot(
            xdata[1:], value[1:] / ydata[1:], color=color, alpha=0.75, label=label
        )

    ax[0].set_ylabel("f(a)")
    ax[0].set_xlabel("a")
    ax[1].set_ylabel("Relative Error")
    ax[1].set_xlabel("a")
    for axis in ax:
        axis.grid()
        axis.legend()

    fig.tight_layout()
    plt.show()


def plot():
    a = np.linspace(-2, 2, 1024)
    y = targetFunc(a)
    plt.plot(a, y)
    plt.grid()
    plt.show()


approx()
