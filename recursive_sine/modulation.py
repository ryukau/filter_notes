"""
Reference:
- https://www.vicanek.de/articles/QuadOsc.pdf
- https://ccrma.stanford.edu/~jos/pasp/Digital_Sinusoid_Generators.html
"""

import matplotlib.pyplot as plt
import numpy as np
import soundfile
from pathlib import Path


def biquad(freqRadian, phaseRadian=0):
    sigS = np.zeros_like(freqRadian)
    sigC = np.zeros_like(freqRadian)

    u0 = np.sin(phaseRadian - freqRadian[0])
    u1 = np.sin(phaseRadian - 2 * freqRadian[0])
    for i in range(len(freqRadian)):
        k = 2 * np.cos(freqRadian[i])
        u2 = u1
        u1 = u0
        u0 = k * u1 - u2

        sigS[i] = u0

    return (sigS, sigC)


def reinsch(freqRadian, phaseRadian=0):
    sigS = np.zeros_like(freqRadian)
    sigC = np.zeros_like(freqRadian)

    u = np.sin(phaseRadian - freqRadian[0])
    v = 2 * np.sin(freqRadian[0] / 2) * np.cos(phaseRadian - freqRadian[0] / 2)
    for i in range(len(freqRadian)):
        A = 2 * np.sin(freqRadian[i] / 2)
        k = A * A

        u = u + v
        v = v - k * u

        sigS[i] = u  # Main output.
        sigC[i] = v

    return (sigS, sigC)


def digitalWaveguide(freqRadian, phaseRadian=0):
    sigS = np.zeros_like(freqRadian)
    sigC = np.zeros_like(freqRadian)

    u = -np.tan(freqRadian[0] / 2) * np.sin(phaseRadian - freqRadian[0])
    v = np.cos(phaseRadian - freqRadian[0])
    for i in range(len(freqRadian)):
        k = np.cos(freqRadian[i])

        s = k * (u + v)
        t = s + u
        u = s - v
        v = t

        sigS[i] = u
        sigC[i] = v  # Main output.

    return (sigS, sigC)


def quadratureStaggered(freqRadian, phaseRadian=0):
    sigS = np.zeros_like(freqRadian)
    sigC = np.zeros_like(freqRadian)

    u = -np.sin(freqRadian[0]) * np.sin(phaseRadian - freqRadian[0])
    v = np.cos(phaseRadian - freqRadian[0])
    for i in range(len(freqRadian)):
        k = np.cos(freqRadian[i])

        t = v
        v = u + k * v
        u = k * v - t

        sigS[i] = u
        sigC[i] = v

    return (sigS, sigC)


def magicCircle(freqRadian, phaseRadian=0):
    sigS = np.zeros_like(freqRadian)
    sigC = np.zeros_like(freqRadian)

    u = np.cos(phaseRadian - freqRadian[0] * 3 / 2)
    v = np.sin(phaseRadian - freqRadian[0])
    for i in range(len(freqRadian)):
        # # Small angle approximation can be used for low frequencies.
        # k = 2 * freqRadian[i]

        k = 2 * np.sin(freqRadian[i] / 2)

        u -= k * v
        v += k * u

        sigS[i] = v
        sigC[i] = u

    return (sigS, sigC)


def coupledForm(freqRadian, phaseRadian=0):
    sigS = np.zeros_like(freqRadian)
    sigC = np.zeros_like(freqRadian)

    u = np.cos(phaseRadian - freqRadian[0])
    v = np.sin(phaseRadian - freqRadian[0])
    for i in range(len(freqRadian)):
        k1 = np.cos(freqRadian[i])
        k2 = np.sin(freqRadian[i])

        u0 = u
        v0 = v
        u = k1 * u0 - k2 * v0
        v = k2 * u0 + k1 * v0

        sigS[i] = v
        sigC[i] = u

    return (sigS, sigC)


def stableQuadrature(freqRadian, phaseRadian=0):
    sigS = np.zeros_like(freqRadian)
    sigC = np.zeros_like(freqRadian)

    u = np.cos(phaseRadian - freqRadian[0])
    v = np.sin(phaseRadian - freqRadian[0])
    for i in range(len(freqRadian)):
        k1 = np.tan(freqRadian[i] / 2)
        k2 = np.sin(freqRadian[i])

        w = u - k1 * v
        v += k2 * w
        u = w - k1 * v

        sigS[i] = v
        sigC[i] = u

    return (sigS, sigC)


def standardMath(freqRadian, phaseRadian=0):
    sigS = np.zeros_like(freqRadian)
    sigC = np.zeros_like(freqRadian)

    freqNormalized = freqRadian / (2 * np.pi)

    phase = 0
    for i in range(len(freqRadian)):
        phase += freqNormalized[i]
        phase -= np.floor(phase)

        sigS[i] = np.sin(2 * np.pi * phase)
        sigC[i] = np.cos(2 * np.pi * phase)

    return (sigS, sigC)


if __name__ == "__main__":
    fs = 48000
    startHz = 10
    peakHz = 20000

    # freqNormalized = np.full(fs, startHz / fs)
    freqNormalized = np.hstack(
        [
            np.geomspace(
                startHz / fs,
                peakHz / fs,
                fs,
            ),
            np.geomspace(
                peakHz / fs,
                startHz / fs,
                fs,
            ),
        ]
    )
    freqRadian = 2 * np.pi * freqNormalized

    def plot(axis, oscFunc):
        sigS, sigC = oscFunc(freqRadian, 0)
        axis.set_title(oscFunc.__name__)
        axis.axhline(1, color="#aaaaaa", lw=0.5, ls="--")
        axis.axhline(-1, color="#aaaaaa", lw=0.5, ls="--")
        axis.plot(sigS, color="red", alpha=0.5, lw=1, label="sin")
        axis.plot(sigC, color="blue", alpha=0.5, lw=1, label="cos")
        axis.grid()
        axis.legend()

        # # Save results as wave files.
        # Path.mkdir("snd", parents=True, exist_ok=True)
        # soundfile.write(f"snd/{oscFunc.__name__}_sin.wav", sigS, fs, subtype="FLOAT")
        # soundfile.write(f"snd/{oscFunc.__name__}_cos.wav", sigC, fs, subtype="FLOAT")

    oscFuncs = [
        biquad,
        reinsch,
        digitalWaveguide,
        quadratureStaggered,
        magicCircle,
        coupledForm,
        stableQuadrature,
        standardMath,
    ]

    nRow = int(np.ceil(len(oscFuncs) / 2))
    fig, ax = plt.subplots(nRow, 2)
    for index, oscFunc in enumerate(oscFuncs):
        col = index // nRow
        row = index % nRow
        plot(ax[row][col], oscFunc)

    fig.suptitle(
        f"Frequency Modulation of Recursive Sine ({startHz} -> {peakHz} -> {startHz} Hz)"
    )
    fig.set_size_inches((10, 1.5 * len(oscFuncs)))
    plt.tight_layout()
    plt.show()

    # Path.mkdir("img", parents=True, exist_ok=True)
    # plt.savefig(f"img/recursive_sine_fm_{startHz}_{peakHz}.svg")
