import sympy
import numpy as np
from scipy.special import lambertw
import matplotlib.pyplot as pyplot
from dataclasses import dataclass


def solveForExpPeakTime():
    a = sympy.Symbol("a", real=True, negative=True)
    d = sympy.Symbol("d", real=True, negative=True)
    t = sympy.Symbol("t", real=True, positive=True)
    E_diff1 = sympy.diff((1 - sympy.exp(a * t)) * sympy.exp(d * t), t)
    solution = sympy.solve(E_diff1, t)
    result = sympy.simplify(solution[0])
    print(sympy.latex(result))


def solvePeakTimeForA():
    a = sympy.Symbol("a", real=True, negative=True)
    d = sympy.Symbol("d", real=True, negative=True)
    t_p = sympy.Symbol("t_p", real=True, positive=True)
    eq = sympy.Eq(t_p, -sympy.log(a / d + 1) / a)
    solution = sympy.solve(eq, a)
    result = sympy.simplify(solution[0])
    print(sympy.latex(result))


def solveForTimeToReachSmall():
    """This couldn't be solved."""
    a = sympy.Symbol("a", real=True, negative=True)
    d = sympy.Symbol("d", real=True, negative=True)
    t = sympy.Symbol("t", real=True, positive=True)
    ε = sympy.Symbol("ε", real=True, positive=True)
    eq = sympy.Eq(ε, (1 - sympy.exp(a * t)) * sympy.exp(d * t))
    solution = sympy.solve(eq, t)
    result = sympy.simplify(solution[0])
    print(sympy.latex(result))


@dataclass
class EnvelopeResult:
    raw: np.array
    normalized: np.array
    peakTime: float


def envelope1(
    attackSeconds: float,
    decaySeconds: float,
    time: np.array,
    eps=np.finfo(np.float64).eps,
):
    """`attackSeconds > 0` and `decaySeconds > 0`."""

    def env(a, d, t):
        return (1 - np.exp(a * t)) * np.exp(d * t)

    a = np.log(eps) / attackSeconds
    d = np.log(eps) / decaySeconds
    peakTime = -np.log1p(a / d) / a
    gain = 1 / env(a, d, peakTime)

    # print(f"E1: t_p={peakTime}, a={a}, d={d}")  # debug

    raw = env(a, d, time)
    return EnvelopeResult(raw, raw * gain, peakTime)


def envelope2(
    attackSeconds: float,
    decaySeconds: float,
    time: np.array,
    eps=np.finfo(np.float64).eps,
):
    """`attackSeconds > 0` and `decaySeconds > 0`."""

    def env(a, d, t):
        return (1 - a**t) * d**t

    a = np.power(eps, 1 / attackSeconds)
    d = np.power(eps, 1 / decaySeconds)
    peakTime = -np.log1p(a / d) / a

    log_a = np.log(a)
    log_d = np.log(d)
    peakTime = np.log(log_d / (log_d + log_a)) / log_a
    gain = 1 / env(a, d, peakTime)

    # print(f"E2: t_p={peakTime}, a={a}, d={d}")  # debug

    raw = env(a, d, time)
    return EnvelopeResult(raw, raw * gain, peakTime)


def envelope3(
    peakSeconds: float,
    releaseSeconds: float,
    time: np.array,
    eps=np.finfo(np.float64).eps,
):
    """
    `peakSeconds > 0`.
    `releaseSeconds` adds to the release time calculated from `peakSeconds`.
    """

    def env(a, d, t):
        return (1 - np.exp(a * t)) * np.exp(d * t)

    D = releaseSeconds - np.log(eps) * peakSeconds
    d = np.log(eps) / D
    x = d * peakSeconds
    a = (lambertw(x * np.exp(x), -1) / peakSeconds - d).real
    gain = 1 / env(a, d, peakSeconds)

    print(f"E3: t_p={peakSeconds}, a={a}, d={d}")  # debug
    print(f"A={-np.log(eps)/np.log(-a)}, D={D}")
    print(np.exp(a))
    print(np.exp(d))

    raw = env(a, d, time)
    return EnvelopeResult(raw, raw * gain, peakSeconds)


def plotEnvelope():
    attackSeconds = 0.1
    decaySeconds = 2.0
    # attackSeconds = 17.06850346467979
    # decaySeconds = 5.604365338911716
    # eps = 1e-5
    eps = np.finfo(np.float64).eps

    samplerate = 48000
    duration = 1
    time = np.linspace(0, duration, duration * samplerate)

    # envFuncs = [envelope1, envelope2, envelope3]
    envFuncs = [envelope3]
    for envFunc in envFuncs:
        result = envFunc(attackSeconds, decaySeconds, time, eps)

        pyplot.figure(figsize=(6, 3))
        pyplot.title(f"ExpAD: A={attackSeconds:3.1f}, D={decaySeconds:3.1f}, ε={eps:e}")
        pyplot.axvline(
            result.peakTime, color="black", alpha=0.5, lw=1, ls=":", label="t_p"
        )
        pyplot.plot(time, result.raw, label="Raw")
        pyplot.plot(time, result.normalized, label="Normalized")
        pyplot.xlabel("Time [s]")
        pyplot.ylabel("Amplitude")
        pyplot.legend()
        pyplot.grid()
        pyplot.tight_layout()
        pyplot.show()


if __name__ == "__main__":
    # solveForExpPeakTime()
    # solvePeakTimeForA()
    # plotTimeToReachSmall()
    plotEnvelope()
