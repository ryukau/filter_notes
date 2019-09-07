import numpy
import matplotlib.pyplot as pyplot

def f(D, a, b, c, e, rho):
    return a * D**3 + b * D**2 + c * D + e - rho

def f_dash(D, a, b, c):
    return 3 * a * D**2 + 2 * b * D + c

def lagrange4(s0, s1, s2, s3):
    return (
        -1 / 6 * s0 + 1 / 2 * s1 - 1 / 2 * s2 + 1 / 6 * s3,
        s0 - 5 / 2 * s1 + 2 * s2 - 1 / 2 * s3,
        -6 / 11 * s0 + 3 * s1 - 3 / 2 * s2 + 1 / 3 * s3,
        s0,
    )

def estimate(sign, rho, s0, s1, s2, s3):
    a, b, c, e = lagrange4(s0 - rho, s1 - rho, s2 - rho, s3 - rho)
    D = 1.5
    for _ in range(128):
        D = D - f(D, a, b, c, e, 0) / f_dash(D, a, b, c)
    return (
        D,
        f(D, a, b, c, e, 0),
        numpy.abs(f_dash(D, a, b, c)),
    )

def blamp_residual(d, slope):
    return (
        slope * (-d**5 / 120 + d**4 / 24 - d**3 / 12 + d**2 / 12 - d / 24 + 1 / 120),
        slope * (d**5 / 40 - d**4 / 12 + d**2 / 3 - d / 2 + 7 / 30),
        slope * (-d**5 / 40 + d**4 / 24 + d**3 / 12 + d**2 / 12 + d / 24 + 1 / 120),
        slope * (d**5 / 120),
    )

def ramp(x, scale):
    return scale * numpy.maximum(0, x)

def plot():
    fix = 1.25
    rho = 1.0
    scale = 1
    sign = 1
    x = numpy.arange(4, dtype=numpy.float64)
    sig = ramp(x - fix, scale) + rho

    x_frac = numpy.linspace(0, len(sig) - 1, 1024)
    true = ramp(x_frac - fix, scale) + rho

    a, b, c, e = lagrange4(*(sig - rho))
    interped = f(x_frac, a, b, c, e, 0) + rho

    D, Dy, slope = estimate(sign, rho, *sig)
    blamp = sig + blamp_residual(D - 1, slope)

    pyplot.title(f"Estimation of $D$ (ρ = {rho})")
    pyplot.plot(x_frac, true, label="Source", color="gray", ls="--", alpha=0.5, zorder=4)
    pyplot.plot(x_frac, interped, label="Lagrange", color="black", zorder=4)
    pyplot.scatter(x, sig, label="Quantized Signal", color="gray", zorder=5)
    pyplot.scatter(x, blamp, label="Result", color="black", marker="x", zorder=5)
    pyplot.scatter(D, rho, label=f"D = {D:.3f}", color="red", marker="+", zorder=5)
    pyplot.scatter(
        fix, rho, label=f"true D = {fix:.3f}", color="blue", marker="^", zorder=5)
    pyplot.grid(zorder=1)
    pyplot.legend()
    pyplot.show()

def plot_edgecase():
    last = 1
    rho = -0.3
    sign = 1
    x = numpy.arange(4, dtype=numpy.float64)
    sig = numpy.array([-0.300, -0.300, -0.001, -0.150])

    x_frac = numpy.linspace(0, len(sig) - 1, 1024)
    a, b, c, e = lagrange4(*(sig - rho))
    interped = f(x_frac, a, b, c, e, 0) + rho

    D, Dy, slope = estimate(sign, rho, *sig)
    blamp = sig + blamp_residual(D - 1, slope)

    pyplot.title(f"Estimation of $D$ (ρ = {rho})")
    pyplot.plot(x_frac, interped, label="Lagrange", color="black", zorder=4)
    pyplot.scatter(x, sig, label="signal", color="gray", zorder=5)
    pyplot.scatter(x, blamp, label="BLAMP", color="black", marker="x", zorder=5)
    pyplot.scatter(D, rho, label=f"D = {D}", color="red", marker="+", zorder=5)
    pyplot.grid(zorder=1)
    pyplot.legend()
    pyplot.show()

# plot()
plot_edgecase()
