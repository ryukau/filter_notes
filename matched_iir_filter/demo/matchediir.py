"""
Notation:
- ω0: center frequency.
- Δω: bandwidth.
- G0: DC gain.
- G1: Nyquist gain.
- G: Peak gain.
- GB: gain at ω ≈ ω0 ± Δω/2.

Reference:
- "Matched Second Order Digital Filters", Martin Vicanek, 14. February 2016.
- "Matched One-Pole Digital Shelving Filters", Martin Vicanek, revised 24. September 2019.

Link to references:
- https://vicanek.de/articles.htm
"""

import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import json
from pathlib import Path


def plotSos(name, targetFreqRadian, sos):
    fs = 2
    worN = np.geomspace(0.0001 * fs, fs / 2, 2048)
    w, h = signal.sosfreqz(sos, worN=worN, fs=fs)

    target = targetFreqRadian * fs / 2 / np.pi
    freq = w
    gain = np.abs(h)
    dB = 20 * np.log10(np.maximum(np.abs(h), 1e-7))
    phase = np.unwrap(np.angle(h))

    fig, ax = plt.subplots(2, 1)
    fig.suptitle(name)

    ax[0].plot(freq, dB)
    ax[0].set_ylabel("Gain [dB]")

    ax[1].plot(freq, phase)
    ax[1].set_ylabel("Phase [rad]")

    for axis in ax:
        axis.grid(which="both", color="#f0f0f0")
        axis.set_xscale("log")
        axis.axvline(target, ls="--", color="black")
    plt.tight_layout()
    plt.plot()
    # plt.show()


def plotComparison(name, targetFreqRadian, sos, impulseResponse):
    fs = 2
    worN = np.geomspace(0.0001 * fs, fs / 2, 2048)
    target = targetFreqRadian * fs / 2 / np.pi

    # Sos.
    w, h = signal.sosfreqz(sos, worN=worN, fs=fs)

    freq = w
    dB = 20 * np.log10(np.maximum(np.abs(h), 1e-7))
    phase = np.unwrap(np.angle(h))

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((8, 6))
    fig.suptitle(name)

    ax[0].plot(freq, dB, color="blue", alpha=0.5, label="python")
    ax[1].plot(freq, phase, color="blue", alpha=0.5, label="python")

    # Impulse response.
    w, h = signal.freqz(impulseResponse, worN=worN, fs=fs)

    freq = w
    dB = 20 * np.log10(np.maximum(np.abs(h), 1e-7))
    phase = np.unwrap(np.angle(h))

    ax[0].plot(freq, dB, color="red", alpha=0.5, label="cpp")
    ax[1].plot(freq, phase, color="red", alpha=0.5, label="cpp")

    for axis in ax:
        axis.set_xscale("log")
        axis.axvline(target, ls="--", color="black", label="f_c")
        axis.grid(which="both", color="#f0f0f0")
        axis.legend()
    ax[0].set_ylabel("Gain [dB]")
    ax[1].set_ylabel("Phase [rad]")
    ax[1].set_xlabel("Frequency [rad/π]")
    plt.tight_layout()
    plt.plot()
    # plt.show()

    imgPath = Path("../img")
    imgPath.mkdir(parents=True, exist_ok=True)
    plt.savefig(imgPath / Path(f"{name}.svg"))
    plt.cla()


def orfanidisPeaking(cutoffRadian, Q, G):
    ω0 = cutoffRadian
    G0 = 1
    G1 = 1
    Δω = ω0 / Q
    GB = np.sqrt(G)

    Ω0 = np.tan(0.5 * ω0)

    G_2 = G * G
    GB_2 = GB * GB
    G0_2 = G0 * G0
    G1_2 = G1 * G1

    W_2 = np.sqrt((G_2 - G1_2) / (G_2 - G0_2)) * Ω0 * Ω0
    ΔΩ = (1 + np.sqrt((GB_2 - G0_2) / (GB_2 - G1_2)) * W_2) * np.tan(0.5 * Δω)

    C = (ΔΩ * ΔΩ) * np.abs(GB_2 - G1_2) - 2 * W_2 * (
        np.abs(GB_2 - G0 * G1) - np.sqrt((GB_2 - G0_2) * (GB_2 - G1_2))
    )
    D = 2 * W_2 * (np.abs(G_2 - G0 * G1) - np.sqrt((G_2 - G0_2) * (G_2 - G1_2)))

    A = np.sqrt((C + D) / np.abs(G_2 - GB_2))
    B = np.sqrt((G_2 * C + GB_2 * D) / np.abs(G_2 - GB_2))

    a0 = 1 + W_2 + A
    a1 = -2 * (1 - W_2) / a0
    a2 = (1 + W_2 - A) / a0
    b0 = (G1 + G0 * W_2 + B) / a0
    b1 = -2 * (G1 - G0 * W_2) / a0
    b2 = (G1 + G0 * W_2 - B) / a0

    return [[b0, b1, b2, 1, a1, a2]]


def massbergLowpass(cutoffRadian, Q):
    ω0 = cutoffRadian

    Q2 = Q * Q
    t1 = np.pi * np.pi / (ω0 * ω0)
    t2 = 1 - t1
    t3 = t1 / Q2
    g1 = 1 / (np.sqrt(t2 * t2 + t3 * t3))

    Ωs = 0
    if Q > np.sqrt(0.5):
        gr = 2 * Q2 / np.sqrt(4 * Q2 - 1)
        ωr = ω0 * np.sqrt(1 - 1 / (2 * Q2))
        Ωs = np.tan(ωr / 2) * np.power((gr * gr - g1 * g1) / (gr * gr - 1), 1 / 4)
    else:
        ωm = ω0 * np.sqrt(
            1 - 1 / 2 / Q2 + np.sqrt((1 - 4 * Q2) / (4 * Q2 * Q2) + 1 / g1)
        )
        Ωs = min(0.5 * ω0 * np.power(1 - g1 * g1, 1 / 4), np.tan(ωm / 2))

    ωz = 2 * np.arctan(Ωs / np.sqrt(g1))
    z_tmp1 = ωz * ωz / (ω0 * ω0)
    z_tmp2 = 1 - z_tmp1
    gz = 1 / (z_tmp2 * z_tmp2 + z_tmp1 / Q2)

    ωp = 2 * np.arctan(Ωs)
    p_tmp1 = ωp * ωp / (ω0 * ω0)
    p_tmp2 = 1 - p_tmp1
    gp = 1 / (p_tmp2 * p_tmp2 + p_tmp1 / Q2)

    gz_2 = gz * gz
    gp_2 = gp * gp
    β = g1 - 1
    Qz = np.sqrt(g1 * g1 * (gp_2 - gz_2) / (gz_2 * (g1 + gp_2) * β * β))
    Qp = np.sqrt(g1 * (gp_2 - gz_2) / ((g1 + gz_2) * β * β))

    Ωs_2 = Ωs * Ωs
    sqrt_g1 = np.sqrt(g1)
    β0 = Ωs_2 + sqrt_g1 * Ωs / Qz + g1
    β1 = 2 * (Ωs_2 - g1)
    β2 = Ωs_2 - sqrt_g1 * Ωs / Qz + g1
    γ = Ωs_2 + Ωs / Qp + 1
    α1 = 2 * (Ωs_2 - 1)
    α2 = Ωs_2 - Ωs / Qp + 1

    b0 = β0 / γ
    b1 = β1 / γ
    b2 = β2 / γ
    a1 = α1 / γ
    a2 = α2 / γ

    return [[b0, b1, b2, 1, a1, a2]]


def solveDenominator(ω0, Q):
    q = 0.5 / Q
    a1 = -2 * np.exp(-q * ω0)
    if q <= 1:
        a1 *= np.cos(np.sqrt(1 - q * q) * ω0)
    else:
        a1 *= np.cosh(np.sqrt(q * q - 1) * ω0)
    a2 = np.exp(-2 * q * ω0)

    sn = np.sin(ω0 / 2)
    φ0 = 1 - sn * sn
    φ1 = sn * sn
    φ2 = 4 * φ0 * φ1

    A0 = (1 + a1 + a2) ** 2
    A1 = (1 - a1 + a2) ** 2
    A2 = -4 * a2

    return (a1, a2, φ0, φ1, φ2, A0, A1, A2)


def matchedLowpass(cutoffRadian, Q):
    ω0 = cutoffRadian

    a1, a2, φ0, φ1, φ2, A0, A1, A2 = solveDenominator(ω0, Q)

    sqrt_B0 = 1 + a1 + a2
    B0 = A0

    R1 = Q * Q * (A0 * φ0 + A1 * φ1 + A2 * φ2)
    B1 = (R1 - B0 * φ0) / φ1

    b0 = 0.5 * (sqrt_B0 + np.sqrt(B1))
    b1 = sqrt_B0 - b0

    return [[b0, b1, 0, 1, a1, a2]]


def matchedHighpass(cutoffRadian, Q):
    ω0 = cutoffRadian

    a1, a2, φ0, φ1, φ2, A0, A1, A2 = solveDenominator(ω0, Q)

    b0 = Q * np.sqrt(A0 * φ0 + A1 * φ1 + A2 * φ2) / (4 * φ1)
    b1 = -2 * b0
    b2 = b0

    return [[b0, b1, b2, 1, a1, a2]]


def matchedBandpass(cutoffRadian, Q):
    ω0 = cutoffRadian

    a1, a2, φ0, φ1, φ2, A0, A1, A2 = solveDenominator(ω0, Q)

    R1 = A0 * φ0 + A1 * φ1 + A2 * φ2
    R2 = -A0 + A1 + 4 * (φ0 - φ1) * A2

    B2 = (R1 - R2 * φ1) / (4 * φ1 * φ1)
    B1 = R2 - 4 * (φ0 - φ1) * B2

    b1 = -0.5 * np.sqrt(B1)
    b0 = 0.5 * (np.sqrt(B2 + b1 * b1) - b1)
    b2 = -b0 - b1

    return [[b0, b1, b2, 1, a1, a2]]


def matchedPeaking(cutoffRadian, Q, gain):
    ω0 = cutoffRadian
    G = gain

    a1, a2, φ0, φ1, φ2, A0, A1, A2 = solveDenominator(ω0, Q)

    R1 = G * G * (A0 * φ0 + A1 * φ1 + A2 * φ2)
    R2 = G * G * (-A0 + A1 + 4 * (φ0 - φ1) * A2)

    B0 = A0
    B2 = (R1 - R2 * φ1 - B0) / (4 * φ1 * φ1)
    B1 = R2 + B0 - 4 * (φ0 - φ1) * B2

    sqrt_B0 = 1 + a1 + a2
    sqrt_B1 = np.sqrt(B1)

    W = 0.5 * (sqrt_B0 + sqrt_B1)
    b0 = 0.5 * (W + np.sqrt(W * W + B2))
    b1 = 0.5 * (sqrt_B0 - sqrt_B1)
    b2 = -B2 / (4 * b0)

    return [[b0, b1, b2, 1, a1, a2]]


def simpleMatchedLowpass(cutoffRadian, Q):
    ω0 = cutoffRadian
    ω0_2 = ω0 * ω0

    a1, a2, _, _, _, _, _, _ = solveDenominator(ω0, Q)

    r0 = 1 + a1 + a2
    r1 = (1 - a1 + a2) * ω0_2 / np.sqrt((1 - ω0_2) ** 2 + ω0_2 / Q / Q)

    b0 = 0.5 * (r0 + r1)
    b1 = r0 - b0

    return [[b0, b1, 0, 1, a1, a2]]


def simpleMatchedHighpass(cutoffRadian, Q):
    ω0 = cutoffRadian
    ω0_2 = ω0 * ω0

    a1, a2, _, _, _, _, _, _ = solveDenominator(ω0, Q)

    r1 = (1 - a1 + a2) / np.sqrt((1 - ω0_2) ** 2 + ω0_2 / Q / Q)

    b0 = 0.25 * r1
    b1 = -2 * b0
    b2 = b0

    return [[b0, b1, b2, 1, a1, a2]]


def simpleMatchedBandpass(cutoffRadian, Q):
    ω0 = cutoffRadian
    ω0_2 = ω0 * ω0

    a1, a2, _, _, _, _, _, _ = solveDenominator(ω0, Q)

    r0 = (1 + a1 + a2) / (ω0 * Q)
    r1 = (1 - a1 + a2) * ω0 / Q / np.sqrt((1 - ω0_2) ** 2 + ω0_2 / Q / Q)

    b0 = 0.5 * r0 + 0.25 * r1
    b1 = -0.5 * r1
    b2 = -b0 - b1

    return [[b0, b1, b2, 1, a1, a2]]


def matchedShelvingOnePole(cutoffRadian, gain):
    fc = cutoffRadian / np.pi
    G = gain

    fm = 0.9
    φm = 1 - np.cos(np.pi * fm)

    def alphabeta(V):
        return 2 / (np.pi * np.pi) * (1 / (fm * fm) + V) - 1 / φm

    α = alphabeta(1 / (G * fc * fc))
    β = alphabeta(G / (fc * fc))

    a1 = -α / (1 + α + np.sqrt(1 + 2 * α))
    b = -β / (1 + β + np.sqrt(1 + 2 * β))
    b0 = (1 + a1) / (1 + b)
    b1 = b * b0

    return [[b0, b1, 0, 1, a1, 0]]


def testResponses(cutoffRadian, qFactor, gainAmp):
    def plot(name, sos):
        plotSos(name, cutoffRadian, sos)

    plot("orfanidisPeaking", orfanidisPeaking(cutoffRadian, qFactor, gainAmp))
    plot("massbergLowpass", massbergLowpass(cutoffRadian, qFactor))
    plot("matchedLowpass", matchedLowpass(cutoffRadian, qFactor))
    plot("matchedHighpass", matchedHighpass(cutoffRadian, qFactor))
    plot("matchedBandpass", matchedBandpass(cutoffRadian, qFactor))
    plot("matchedPeaking", matchedPeaking(cutoffRadian, qFactor, gainAmp))
    plot("simpleMatchedLowpass", simpleMatchedLowpass(cutoffRadian, qFactor))
    plot("simpleMatchedHighpass", simpleMatchedHighpass(cutoffRadian, qFactor))
    plot("simpleMatchedBandpass", simpleMatchedBandpass(cutoffRadian, qFactor))
    plot("matchedShelvingOnePole", matchedShelvingOnePole(cutoffRadian, gainAmp))

    plt.show()


def compareResponsesToCpp():
    with open("impulse_response.json", "r", encoding="utf-8") as fi:
        data = json.load(fi)

    cutoff = 2 * np.pi * data["cutoffNormalized"]  # Cutoff frquency [rad].
    reso = data["Q"]  # Quality factor, or resonance.
    gainAmp = 10 ** (data["gainDecibel"] / 20)  # Gain in amplitude.

    func = {
        "orfanidisPeaking": lambda fc, Q, G: orfanidisPeaking(fc, Q, G),
        "massbergLowpass": lambda fc, Q, G: massbergLowpass(fc, Q),
        "matchedLowpass": lambda fc, Q, G: matchedLowpass(fc, Q),
        "matchedHighpass": lambda fc, Q, G: matchedHighpass(fc, Q),
        "matchedBandpass": lambda fc, Q, G: matchedBandpass(fc, Q),
        "matchedPeaking": lambda fc, Q, G: matchedPeaking(fc, Q, G),
        "simpleMatchedLowpass": lambda fc, Q, G: simpleMatchedLowpass(fc, Q),
        "simpleMatchedHighpass": lambda fc, Q, G: simpleMatchedHighpass(fc, Q),
        "simpleMatchedBandpass": lambda fc, Q, G: simpleMatchedBandpass(fc, Q),
        "matchedShelvingOnePole": lambda fc, Q, G: matchedShelvingOnePole(fc, G),
    }

    irs = data["impulseResponse"]
    for key, ir in irs.items():
        plotComparison(key, cutoff, func[key](cutoff, reso, gainAmp), ir)

    # # For debugging single filter.
    # key = "massbergLowpass"
    # plotComparison(key, cutoff, func[key](cutoff, reso, gainAmp), irs[key])

    # plt.show()


if __name__ == "__main__":
    # testResponses(np.pi * 0.01, np.sqrt(2) / 2, 10)
    compareResponsesToCpp()
