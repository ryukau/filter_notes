import numpy as np
import matplotlib.pyplot as plt


def biquad(freqNormalized, initialPhase):
    """
    `freqNormalized` is in [0, 0.5).
    `initialPhase` is in [0, 1).

    Return value is `(sin, cos)`.
    """
    omega = 2 * np.pi * freqNormalized
    phi = 2 * np.pi * initialPhase
    u1 = np.sin(phi - omega)
    u2 = np.sin(phi - 2 * omega)
    k = 2 * np.cos(omega)
    while True:
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0
        yield (u0, 0)


def reinsch(freqNormalized, initialPhase):
    omega = 2 * np.pi * freqNormalized
    phi = 2 * np.pi * initialPhase
    A = 2 * np.sin(omega / 2)
    u = np.sin(phi - omega)
    v = A * np.cos(phi - omega / 2)
    k = A * A
    while True:
        u = u + v
        v = v - k * u
        yield (u, v)  # `u` is main output.


def digitalWaveguide(freqNormalized, initialPhase):
    omega = 2 * np.pi * freqNormalized
    phi = 2 * np.pi * initialPhase
    u = -np.tan(omega / 2) * np.sin(phi - omega)
    v = np.cos(phi - omega)
    k = np.cos(omega)
    while True:
        s = k * (u + v)
        t = s + u
        u = s - v
        v = t
        yield (u, v)  # `v` is main output.


def quadratureStaggered(freqNormalized, initialPhase):
    omega = 2 * np.pi * freqNormalized
    phi = 2 * np.pi * initialPhase
    u = -np.sin(omega) * np.sin(phi - omega)
    v = np.cos(phi - omega)
    k = np.cos(omega)
    while True:
        t = v
        v = u + k * v
        u = k * v - t
        yield (u, v)  # `v` is main output.


def magicCircle(freqNormalized, initialPhase):
    omega = 2 * np.pi * freqNormalized
    phi = 2 * np.pi * initialPhase
    u = np.cos(phi - omega * 3 / 2)
    v = np.sin(phi - omega)
    k = 2 * np.sin(omega / 2)
    while True:
        u -= k * v
        v += k * u
        yield (v, u)


def coupledForm(freqNormalized, initialPhase):
    omega = 2 * np.pi * freqNormalized
    phi = 2 * np.pi * initialPhase
    u = np.cos(phi - omega)
    v = np.sin(phi - omega)
    k1 = np.cos(omega)
    k2 = np.sin(omega)
    while True:
        u0 = u
        v0 = v
        u = k1 * u0 - k2 * v0
        v = k2 * u0 + k1 * v0
        yield (v, u)


def stableQuadrature(freqNormalized, initialPhase):
    omega = 2 * np.pi * freqNormalized
    phi = 2 * np.pi * initialPhase
    u = np.cos(phi - omega)
    v = np.sin(phi - omega)
    k1 = np.tan(omega / 2)
    k2 = np.sin(omega)
    while True:
        w = u - k1 * v
        v += k2 * w
        u = w - k1 * v
        yield (v, u)


def plotInitialPhase(oscFunc):
    sampleRate = 1024
    freqHz = 1
    nPhase = 6

    cmap = plt.get_cmap("magma")
    for index, initialPhase in enumerate(np.linspace(0, 1, nPhase, endpoint=False)):
        oscGen = oscFunc(freqHz / sampleRate, initialPhase)
        sig = np.zeros(1024)
        for i in range(len(sig)):
            sig[i] = next(oscGen)[0]
        plt.plot(sig, alpha=0.5, color=cmap(index / nPhase), label=f"{initialPhase:.3}")

    plt.grid()
    plt.legend()
    plt.show()


def plotSincos(oscFunc, initialPhase=0, normalize=False):
    sampleRate = 1024
    freqHz = 1

    oscGen = oscFunc(freqHz / sampleRate, initialPhase)

    length = 1 * sampleRate
    sigS = np.zeros(length)
    sigC = np.zeros(length)
    for i in range(length):
        output = next(oscGen)
        sigS[i] = output[0]
        sigC[i] = output[1]

    def normalizeAmp(sig):
        max = np.max(np.abs(sig))
        if max == 0:
            return sig
        return sig / max

    if normalize:
        sigS = normalizeAmp(sigS)
        sigC = normalizeAmp(sigC)

    # # Print first value.
    # refS = np.sin(2 * np.pi * initialPhase)
    # refC = np.cos(2 * np.pi * initialPhase)
    # print(f"--- {oscFunc.__name__}")
    # print(f"output (sin      , cos      ) : ({sigS[0]}, {sigC[0]}")
    # print(f"diff   (sin - ref, cos - ref) : ({sigS[0] - refS}, {sigC[0] - refC}")

    plt.title(oscFunc.__name__)
    plt.plot(sigS, alpha=0.5, color="red", label=f"Sin")
    plt.plot(sigC, alpha=0.5, color="blue", label=f"Cos")

    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    initialPhase = 0
    plotSincos(biquad, initialPhase, False)
    plotSincos(reinsch, initialPhase, False)
    plotSincos(digitalWaveguide, initialPhase, False)
    plotSincos(quadratureStaggered, initialPhase, False)
    plotSincos(magicCircle, initialPhase, False)
    plotSincos(coupledForm, initialPhase, False)
    plotSincos(stableQuadrature, initialPhase, False)
