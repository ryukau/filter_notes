"""
References:

- Parker, J. D., Zavalishin, V., & Le Bivic, E. (2016, September). Reducing the aliasing of nonlinear waveshaping using continuous-time convolution. In Proc. Int. Conf. Digital Audio Effects (DAFx-16), Brno, Czech Republic (pp. 137-144).
- Bilbao, S., Esqueda, F., Parker, J. D., & Välimäki, V. (2017). Antiderivative antialiasing for memoryless nonlinearities. IEEE Signal Processing Letters, 24(7), 1049-1053.
- La Pastina, P. P., D'Angelo, S., & Gabrielli, L. (2021, September). Arbitrary-order IIR antiderivative antialiasing. In 2021 24th International Conference on Digital Audio Effects (DAFx) (pp. 9-16). IEEE.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.special as special
import soundfile


def generateSine(durationSamples: int, frequencyNormalized: float) -> np.ndarray:
    phase = 0
    sig = np.empty(durationSamples)
    for i in range(len(sig)):
        sig[i] = np.sin(2 * np.pi * phase)
        phase += frequencyNormalized
        phase -= np.floor(phase)
    return sig


def hardclipJ0(x, gain=1):
    return np.clip(gain * x, -1, 1)


def hardclipJ1(x):
    absed = np.abs(x)
    if absed < 1:
        return x * x / 2
    return absed - 0.5


def hardclipJ2(x):
    if np.abs(x) < 1:
        return x * x * x / 6
    return (x * x / 2 + 1 / 6) * np.sign(x) - (x / 2)


def hardclipJ3(x):
    if np.abs(x) < 1:
        return x * x * x * x / 24
    return (x * x * x / 6 + x / 6) * np.sign(x) - (x * x / 4)


def hardclipParkerJ2(x):
    if np.abs(x) < 1:
        return x * x * x / 3
    return ((x * x / 2) + 1 / 6) * np.sign(x)


def tanhParkerJ1(x):
    return np.log(np.cosh(x))


def tanhParkerJ2(x):
    alpha = np.exp(-2 * x)
    return 0.5 * (x * (x + 2 * np.log(alpha + 1)) - special.spence(-alpha))


def adaaParker2(x, f0, f1, f2):
    tolerance = np.finfo(np.float64).eps
    y = np.zeros_like(x)

    x1 = 0
    x2 = 0
    p1 = 0
    p2 = 0
    q1 = 0
    q2 = 0
    for n, _ in enumerate(x):
        x0 = x[n]
        p0 = f1(x0)
        q0 = f2(x0)

        d0 = x0 - x1
        t0 = 0
        if np.abs(d0) < tolerance:
            t0 = 0.5 * f0((x0 + 2 * x1) / 3)
        else:
            t0 = (x0 * (p0 - p1) - (q0 - q1)) / (d0 * d0)

        d1 = x2 - x1
        t1 = 0
        if np.abs(d1) < tolerance:
            t1 = 0.5 * f0((x2 + 2 * x1) / 3)
        else:
            t1 = (x2 * (p2 - p1) - (q2 - q1)) / (d1 * d1)

        y[n] = t0 + t1

        p2 = p1
        p1 = p0
        q2 = q1
        q1 = q0
        x2 = x1
        x1 = x0
    return y


def adaaBilbao1(x, f0, f1):
    tolerance = 1 / 2**24
    y = np.zeros_like(x)
    x1 = 0
    for n, _ in enumerate(x):
        if np.abs(x[n] - x1) < tolerance:  # fallback
            y[n] = f0((x[n] + x1) / 2)
        else:
            y[n] = (f1(x[n]) - f1(x1)) / (x[n] - x1)
        x1 = x[n]
    return y


def adaaBilbao2(x, f0, f1, f2):
    tolerance = 1 / 2**16
    y = np.zeros_like(x)

    def calcD(x0, x1):
        if np.abs(x0 - x1) < tolerance:
            return f1((x0 + x1) / 2)
        return (f2(x0) - f2(x1)) / (x0 - x1)

    x1 = 0
    x2 = 0
    s1 = 0
    for n, _ in enumerate(x):
        x0 = x[n]
        s0 = calcD(x0, x1)
        if np.abs(x0 - x2) < tolerance:  # fallback
            x_bar = (x0 + x2) / 2
            delta = x_bar - x1
            if delta < tolerance:
                y[n] = f0((x_bar + x0) / 2)
            else:
                y[n] = (2 / delta) * (f1(x_bar) + (f2(x1) - f2(x_bar)) / delta)
            pass
        else:
            y[n] = (2 / (x0 - x2)) * (calcD(x0, x1) - s1)
        s1 = s0
        x2 = x1
        x1 = x0
    return y


# computing the integral
def integrateAaiir(x_z1, x, s):
    if x_z1 < x:
        if x <= -1:
            return (1 - np.exp(s)) / s
        elif x > -1 and x <= 1 and x_z1 <= -1:
            return (
                (x - x_z1) * np.exp(s * (x + 1) / (x - x_z1))
                - (s + 1) * x
                + x_z1
                - s * np.exp(s)
            ) / (s * s)
        elif x > 1 and x_z1 <= -1:
            return (
                (x - x_z1)
                * (np.exp(s * (x + 1) / (x - x_z1)) - np.exp(s * (x - 1) / (x - x_z1)))
                - s * (np.exp(s) + 1)
            ) / (s * s)
        elif x_z1 > -1 and x <= 1:
            return (np.exp(s) * (x + (s - 1) * x_z1) - (s + 1) * x + x_z1) / (s * s)
        elif x_z1 > -1 and x_z1 <= 1 and x > 1:
            return (
                (x_z1 - x) * np.exp(s * (x - 1) / (x - x_z1))
                + np.exp(s) * (x + (s - 1) * x_z1)
                - s
            ) / (s * s)
        else:
            return (np.exp(s) - 1) / s
    else:
        if x_z1 <= -1:
            return (1 - np.exp(s)) / s
        elif x_z1 > -1 and x_z1 <= 1 and x <= -1:
            return (
                (x_z1 - x) * np.exp(s * (x + 1) / (x - x_z1))
                + np.exp(s) * (x + (s - 1) * x_z1)
                + s
            ) / (s * s)
        elif x_z1 > 1 and x <= -1:
            return (
                (x_z1 - x)
                * (
                    np.exp((s * (x + 1)) / (x - x_z1))
                    - np.exp((s * (x - 1)) / (x - x_z1))
                )
                + s * (np.exp(s) + 1)
            ) / (s * s)
        elif x > -1 and x_z1 <= 1:
            return (np.exp(s) * (x + (s - 1) * x_z1) - (s + 1) * x + x_z1) / (s * s)
        elif x > -1 and x <= 1 and x_z1 > 1:
            return (
                (x - x_z1) * np.exp((s * (x - 1)) / (x - x_z1))
                - (s + 1) * x
                + x_z1
                + s * np.exp(s)
            ) / (s * s)
        else:
            return (np.exp(s) - 1) / s


def aaiirHardComplex(x, B, beta, tol):
    # precalculating some values
    E = np.exp(beta)
    F = (E - 1) / beta

    # initial conditions
    x_z1 = 0
    Y_z1 = 0
    y = np.zeros_like(x)

    # audio processing loop
    for n in range(len(x)):
        if np.abs(x[n] - x_z1) < tol:
            I = hardclipJ0(0.5 * (x[n] + x_z1)) * F
        else:
            I = integrateAaiir(x_z1, x[n], beta)
        Y = E * Y_z1 + 2 * B * I
        y[n] = Y.real

        x_z1 = x[n]
        Y_z1 = Y
    return y


def hardclipAaiir(x, b, a, r, p):
    y = np.zeros_like(x)
    for ord in range(0, len(r), 2):
        # poles are in conjg pairs, must take only one for each
        ri = r[ord]
        zi = p[ord]
        y = y + aaiirHardComplex(x, ri, zi, 1e-9)
    return y


def hardclipParker2(x, gain=10):
    return adaaParker2(gain * x, hardclipJ0, hardclipJ1, hardclipParkerJ2)


def hardclipBilbao1(x, gain=10):
    return adaaBilbao1(gain * x, hardclipJ0, hardclipJ1)


def hardclipBilbao2(x, gain=10):
    return adaaBilbao2(gain * x, hardclipJ0, hardclipJ1, hardclipJ2)


def hardclipAaiirButter(x, gain=10):
    order = 4
    cutoffRadian = 2 * np.pi * 18000 / 48000

    b, a = signal.butter(order, cutoffRadian, "low", analog=True)
    r, p, k = signal.residue(b, a)

    return hardclipAaiir(gain * x, b, a, r, p)


def hardclipAaiirChebyII(x, gain=10):
    order = 8
    stopbandDecibel = 60
    cutoffRadian = 2 * np.pi * 22000 / 48000

    b, a = signal.cheby2(order, stopbandDecibel, cutoffRadian, "low", analog=True)
    r, p, k = signal.residue(b, a)

    return hardclipAaiir(gain * x, b, a, r, p)


if __name__ == "__main__":
    sampleRate = 48000

    frequencyHz = 61
    durationSecond = 1
    gain = 100
    clipFunctions = [
        hardclipJ0,
        hardclipParker2,
        hardclipAaiirButter,
        hardclipAaiirChebyII,
        hardclipBilbao1,
        hardclipBilbao2,
    ]

    for clipFn in clipFunctions:
        sig = generateSine(int(sampleRate * durationSecond), frequencyHz / sampleRate)
        # sig *= signal.windows.tukey(len(sig), 0.1 / durationSecond)
        sig = clipFn(sig, gain)

        spc = np.fft.rfft(sig, norm="forward")
        absed = np.abs(spc)
        mag = 20 * np.log10(absed)
        freq = np.fft.rfftfreq(len(sig), 1 / sampleRate)

        soundfile.write(f"{clipFn.__name__}.wav", 0.1 * sig, sampleRate, "FLOAT")

        fig, ax = plt.subplots(2, 1)
        fig.set_size_inches(8, 6)
        fig.suptitle(f"{clipFn.__name__}")

        ax[0].plot(sig)
        ax[0].set_xlabel("Time [sample]")
        ax[0].set_ylabel("Amplitude")

        ax[1].plot(freq, mag)
        # ax[1].set_xscale("log")
        ax[1].set_xlabel("Frequency [Hz]")
        ax[1].set_ylabel("Gain [dB]")
        ax[1].set_xlim([10, sampleRate // 2])
        ax[1].set_ylim([-120, 10])

        for axis in ax:
            axis.grid()
        plt.tight_layout()
        plt.show()
