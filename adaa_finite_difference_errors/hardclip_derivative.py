import numpy as np
import matplotlib.pyplot as plt
import math


def f0(x):
    if x <= -1.0:
        return -1.0
    if x >= 1.0:
        return 1.0
    return x


def f1(x):
    z = abs(x)
    if z < 1.0:
        return 0.5 * x * x
    return z - 0.5


def f2(x):
    if abs(x) < 1.0:
        return (x * x * x) / 6.0
    s = 1.0 if x >= 0 else -1.0
    return (s * x - 1) * x / 2 + s / 6


def f3(x):
    x2 = x * x
    if x2 < 1.0:
        return (x2 * x2) / 24.0
    s = 1.0 if x >= 0 else -1.0
    return (((s * 4 * x - 6) * x + s * 4) * x - 1) / 24


def f4(x):
    z = abs(x)
    if z < 1.0:
        return (x**5) / 120.0
    s = 1.0 if x >= 0 else -1.0
    return s * (((((z - 2.0) * z + 2.0) * z - 1.0) * z / 24.0) + 1.0 / 120.0)


def g1(x):
    if abs(x) >= 1.0:
        return 0.0
    return 1.0


def g2(x):
    return 0.0


def plot_results():
    x = np.linspace(-8, 8, 2048)

    fig, axes = plt.subplots(7, 1, figsize=(8, 8), sharex=True)

    plots = [
        ([g2(v) for v in x], "g2", None),
        ([g1(v) for v in x], "g1", None),
        ([f0(v) for v in x], "f0", None),
        ([f1(v) for v in x], "f1", None),
        ([f2(v) for v in x], "f2", None),
        ([f3(v) for v in x], "f3", None),
        ([f4(v) for v in x], "f4", None),
    ]

    for ax, (y, title, ref_func) in zip(axes, plots):
        ax.plot(x, y, "b-", lw=1.5, label="Exact")
        if ref_func:
            ax.plot(x, ref_func(x), "g--", alpha=0.5, label="Reference")
            ax.legend(loc="upper left")
        ax.set_title(title)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    # plt.show()


def plot_diff():
    """
    offset of `g1` changes with the starting point of `x`.
    """
    x = np.linspace(-8, 8, 2048)
    dx = x[1] - x[0]

    g1_y = np.cumsum([g1(v) for v in x]) * dx

    plt.figure()
    plt.plot(x, g1_y - np.mean(g1_y), label="g1")
    plt.plot(x, [f0(v) for v in x], label="f0")
    plt.plot(x[:-1], np.diff([f1(v) for v in x], 1) / dx**1, label="f1")
    plt.plot(x[1:-1], np.diff([f2(v) for v in x], 2) / dx**2, label="f2")
    plt.plot(x[1:-2], np.diff([f3(v) for v in x], 3) / dx**3, label="f3")
    plt.plot(x[2:-2], np.diff([f4(v) for v in x], 4) / dx**4, label="f4")

    plt.legend()
    plt.grid()
    # plt.show()


def verify_smoothness():
    x = np.linspace(9.0, 11.0, 2000)
    dx = x[1] - x[0]
    d4 = np.diff([f4(v) for v in x], 4) / dx**4
    x_valid = x[2:-2]

    plt.figure(figsize=(10, 6))
    plt.plot(x_valid, d4, "r-", linewidth=2, label="Recovered f0 (Numerical)")
    plt.plot(x, [f0(v) for v in x], "k--", alpha=0.5, label="Ideal f0")
    plt.title("Robust Verification (Edges Trimmed)")
    plt.xlabel("x")
    plt.legend()
    plt.grid(True, alpha=0.3)

    # plt.show()


if __name__ == "__main__":
    plot_results()
    plot_diff()
    verify_smoothness()
    plt.show()
