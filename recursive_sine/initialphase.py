import numpy as np
import matplotlib.pyplot as plt


def biquad(freqNormalized, initialPhase):
    """
    `freqNormalized` is in [0, 0.5).
    `initialPhase` is in [0, 1).
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
        yield u0


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
        yield u


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
        yield v


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
        yield v


def magicCircle(freqNormalized, initialPhase):
    omega = 2 * np.pi * freqNormalized
    phi = 2 * np.pi * initialPhase
    u = np.cos(phi - omega * 3 / 2)
    v = np.sin(phi - omega)
    k = 2 * np.sin(omega / 2)
    while True:
        u -= k * v
        v += k * u
        yield v


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
        yield v


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
        yield v


sampleRate = 1024
freqHz = 1
nPhase = 6

cmap = plt.get_cmap("magma")
for index, initialPhase in enumerate(np.linspace(0, 1, nPhase, endpoint=False)):
    bq = stableQuadrature(freqHz / sampleRate, initialPhase)
    sig = np.empty(1024)
    for i in range(len(sig)):
        sig[i] = next(bq)
    plt.plot(sig, alpha=0.5, color=cmap(index / nPhase), label=f"{initialPhase:.3}")

plt.grid()
plt.legend()
plt.show()
