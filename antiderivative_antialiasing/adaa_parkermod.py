"""
A reference implemetation of antiderivative antialiasing using cosine window in continuous domain. The expressions were not simplified.
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


def normalize(sig):
    peak = np.max(np.abs(sig))
    if peak <= np.finfo(np.float64).eps:
        return sig
    return sig / peak


def halfrectJ0(x, gain=1):
    return np.maximum(0, gain * x)


def halfrectJ1(x):
    if x < 0:
        return 0
    return x * x / 2


def halfrectJ2(x):
    if x < 0:
        return x * x * x / 3


def hardclipJ0(x, gain=1):
    return np.clip(gain * x, -1, 1)


def hardclipJ1(x):
    absed = np.abs(x)
    if np.abs(x) < 1:
        return x * x / 2
    return absed - 0.5


def hardclipJ2(x):
    if np.abs(x) < 1:
        return x * x * x / 6
    return ((x * x / 2) + (1 / 6)) * np.sign(x) - (x / 2)


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
    tolerance = 1 / 2**24
    y = np.zeros_like(x)

    x1 = 0
    x2 = 0
    for n, _ in enumerate(x):
        x0 = x[n]

        d0 = x0 - x1
        t0 = 0
        if np.abs(d0) < tolerance:
            t0 = 0.5 * f0((x0 + 2 * x1) / 3)
        else:
            t0 = (x0 * (f1(x0) - f1(x1)) - (f2(x0) - f2(x1))) / (d0 * d0)

        d1 = x1 - x2
        t1 = 0
        if np.abs(d1) < tolerance:
            t1 = 0.5 * f0((x2 + 2 * x1) / 3)
        else:
            t1 = (x2 * (f1(x2) - f1(x1)) - (f2(x2) - f2(x1))) / (d1 * d1)

        y[n] = t0 + t1

        x2 = x1
        x1 = x0
    return y


def adaaCosineHalfRectJ2(x, gain=1):
    tolerance = 1 / 2**24
    x = gain * x.copy()
    y = np.zeros_like(x)

    x1 = 0
    x2 = 0
    for n, _ in enumerate(x):
        x0 = x[n]

        d0 = x0 - x1
        t0 = 0
        if np.abs(d0) < tolerance:
            mid = 0.5 + 2 / (np.pi * np.pi)
            t0 = 0.5 * halfrectJ0((x2 + mid * (x1 - x2)))
        else:
            if x0 < 0:
                if x1 < 0:
                    t0 = 0
                else:
                    t0 = (
                        -(x0**2) * np.cos(np.pi * x0 / (x0 - x1)) / 2
                        - x0**2 / 2
                        + x0 * x1 * np.cos(np.pi * x0 / (x0 - x1))
                        + x0 * x1
                        - x1**2 * np.cos(np.pi * x0 / (x0 - x1)) / 2
                        - np.pi**2 * x1**2 / 4
                        - x1**2 / 2
                    ) / (np.pi**2 * (x0 - x1))
            else:
                if x1 < 0:
                    t0 = (
                        np.pi**2 * x0**2 * (-x0 + x1)
                        + 2 * np.pi**2 * x0**2 * (x0 - x1)
                        + 2
                        * (x0 - x1) ** 2
                        * (
                            x0 * np.cos(np.pi * x0 / (x0 - x1))
                            - x0
                            - x1 * np.cos(np.pi * x0 / (x0 - x1))
                            + x1
                        )
                    ) / (4 * np.pi**2 * (x0 - x1) ** 2)
                else:
                    t0 = (-x0 + x1 + np.pi**2 * (x0 + x1) / 4) / np.pi**2

        d1 = x1 - x2
        t1 = 0
        if np.abs(d1) < tolerance:
            mid = 0.5 - 2 / (np.pi * np.pi)
            t1 = 0.5 * halfrectJ0((x1 + mid * (x2 - x1)))
        else:
            if x1 < 0:
                if x2 < 0:
                    t1 = 0
                else:
                    t1 = (
                        -np.pi**2 * x1**2
                        + np.pi**2 * (x1 - x2) * (x1 + x2)
                        + 2
                        * (x1 - x2)
                        * (
                            x1 * np.cos(np.pi * x1 / (x1 - x2))
                            + x1
                            - x2 * np.cos(np.pi * x1 / (x1 - x2))
                            - x2
                        )
                    ) / (4 * np.pi**2 * (x1 - x2))
            else:
                if x2 < 0:
                    t1 = (
                        -(x1**2) * np.cos(np.pi * x1 / (x1 - x2)) / 2
                        + x1**2 / 2
                        + np.pi**2 * x1**2 / 4
                        + x1 * x2 * np.cos(np.pi * x1 / (x1 - x2))
                        - x1 * x2
                        - x2**2 * np.cos(np.pi * x1 / (x1 - x2)) / 2
                        + x2**2 / 2
                    ) / (np.pi**2 * (x1 - x2))
                else:
                    t1 = (x1 - x2 + np.pi**2 * (x1 + x2) / 4) / np.pi**2

        y[n] = t0 + t1

        x2 = x1
        x1 = x0
    return y


def adaaCosineHardclipJ2(x, gain=1):
    tolerance = 1 / 2**24
    x = gain * x.copy()
    y = np.zeros_like(x)

    x1 = 0
    x2 = 0
    for n, _ in enumerate(x):
        x0 = x[n]

        d0 = x0 - x1
        t0 = 0
        if np.abs(d0) < tolerance:
            mid = 0.5 + 2 / (np.pi * np.pi)
            t0 = 0.5 * hardclipJ0((x0 + mid * (x1 - x0)))
        else:
            if x0 < -1 and x1 < -1:
                t0 = -1 / 2
            elif x0 < -1 and abs(x1) <= 1:
                t0 = (
                    -2 * x0**2 * np.cos(np.pi * (x0 + 1) / (x0 - x1))
                    - 2 * x0**2
                    + 4 * x0 * x1 * np.cos(np.pi * (x0 + 1) / (x0 - x1))
                    + 4 * x0 * x1
                    - 2 * np.pi**2 * x0
                    - 2 * x1**2 * np.cos(np.pi * (x0 + 1) / (x0 - x1))
                    - np.pi**2 * x1**2
                    - 2 * x1**2
                    - np.pi**2
                ) / (4 * np.pi**2 * x0 - 4 * np.pi**2 * x1)
            elif x0 < -1 and x1 > 1:
                t0 = (
                    2
                    * x0**2
                    * np.sin(np.pi / (x0 - x1))
                    * np.sin(np.pi * x0 / (x0 - x1))
                    - 4
                    * x0
                    * x1
                    * np.sin(np.pi / (x0 - x1))
                    * np.sin(np.pi * x0 / (x0 - x1))
                    - np.pi**2 * x0
                    + 2
                    * x1**2
                    * np.sin(np.pi / (x0 - x1))
                    * np.sin(np.pi * x0 / (x0 - x1))
                    - np.pi**2 * x1
                ) / (2 * np.pi**2 * x0 - 2 * np.pi**2 * x1)
            elif abs(x0) <= 1 and x1 < -1:
                t0 = (
                    2 * x0**2 * np.cos(np.pi * (x0 + 1) / (x0 - x1))
                    - 2 * x0**2
                    + np.pi**2 * x0**2
                    - 4 * x0 * x1 * np.cos(np.pi * (x0 + 1) / (x0 - x1))
                    + 4 * x0 * x1
                    + 2 * x1**2 * np.cos(np.pi * (x0 + 1) / (x0 - x1))
                    - 2 * x1**2
                    + 2 * np.pi**2 * x1
                    + np.pi**2
                ) / (4 * np.pi**2 * x0 - 4 * np.pi**2 * x1)
            elif abs(x0) <= 1 and abs(x1) <= 1:
                t0 = (-4 * x0 + np.pi**2 * x0 + 4 * x1 + np.pi**2 * x1) / (4 * np.pi**2)
            elif abs(x0) <= 1 and x1 > 1:
                t0 = (
                    2 * x0**2 * np.cos(np.pi * (x0 - 1) / (x0 - x1))
                    - 2 * x0**2
                    + np.pi**2 * x0**2
                    - 4 * x0 * x1 * np.cos(np.pi * (x0 - 1) / (x0 - x1))
                    + 4 * x0 * x1
                    + 2 * x1**2 * np.cos(np.pi * (x0 - 1) / (x0 - x1))
                    - 2 * x1**2
                    - 2 * np.pi**2 * x1
                    + np.pi**2
                ) / (4 * np.pi**2 * x0 - 4 * np.pi**2 * x1)
            elif x0 > 1 and x1 < -1:
                t0 = (
                    -2
                    * x0**2
                    * np.sin(np.pi / (x0 - x1))
                    * np.sin(np.pi * x0 / (x0 - x1))
                    + 4
                    * x0
                    * x1
                    * np.sin(np.pi / (x0 - x1))
                    * np.sin(np.pi * x0 / (x0 - x1))
                    + np.pi**2 * x0
                    - 2
                    * x1**2
                    * np.sin(np.pi / (x0 - x1))
                    * np.sin(np.pi * x0 / (x0 - x1))
                    + np.pi**2 * x1
                ) / (2 * np.pi**2 * x0 - 2 * np.pi**2 * x1)
            elif x0 > 1 and abs(x1) <= 1:
                t0 = (
                    -2 * x0**2 * np.cos(np.pi * (x0 - 1) / (x0 - x1))
                    - 2 * x0**2
                    + 4 * x0 * x1 * np.cos(np.pi * (x0 - 1) / (x0 - x1))
                    + 4 * x0 * x1
                    + 2 * np.pi**2 * x0
                    - 2 * x1**2 * np.cos(np.pi * (x0 - 1) / (x0 - x1))
                    - np.pi**2 * x1**2
                    - 2 * x1**2
                    - np.pi**2
                ) / (4 * np.pi**2 * x0 - 4 * np.pi**2 * x1)
            elif x0 > 1 and x1 > 1:
                t0 = 1 / 2

        d1 = x1 - x2
        t1 = 0
        if np.abs(d1) < tolerance:
            mid = 0.5 - 2 / (np.pi * np.pi)
            t1 = 0.5 * hardclipJ0((x1 + mid * (x2 - x1)))
        else:
            if x1 < -1 and x2 < -1:
                t1 = -1 / 2
            elif x1 < -1 and abs(x2) <= 1:
                t1 = (
                    2 * x1**2 * np.cos(np.pi * (x1 + 1) / (x1 - x2))
                    + 2 * x1**2
                    - 4 * x1 * x2 * np.cos(np.pi * (x1 + 1) / (x1 - x2))
                    - 4 * x1 * x2
                    - 2 * np.pi**2 * x1
                    + 2 * x2**2 * np.cos(np.pi * (x1 + 1) / (x1 - x2))
                    - np.pi**2 * x2**2
                    + 2 * x2**2
                    - np.pi**2
                ) / (4 * np.pi**2 * x1 - 4 * np.pi**2 * x2)
            elif x1 < -1 and x2 > 1:
                t1 = (
                    -2
                    * x1**2
                    * np.sin(np.pi / (x1 - x2))
                    * np.sin(np.pi * x1 / (x1 - x2))
                    + 4
                    * x1
                    * x2
                    * np.sin(np.pi / (x1 - x2))
                    * np.sin(np.pi * x1 / (x1 - x2))
                    - np.pi**2 * x1
                    - 2
                    * x2**2
                    * np.sin(np.pi / (x1 - x2))
                    * np.sin(np.pi * x1 / (x1 - x2))
                    - np.pi**2 * x2
                ) / (2 * np.pi**2 * x1 - 2 * np.pi**2 * x2)
            elif abs(x1) <= 1 and x2 < -1:
                t1 = (
                    -2 * x1**2 * np.cos(np.pi * (x1 + 1) / (x1 - x2))
                    + 2 * x1**2
                    + np.pi**2 * x1**2
                    + 4 * x1 * x2 * np.cos(np.pi * (x1 + 1) / (x1 - x2))
                    - 4 * x1 * x2
                    - 2 * x2**2 * np.cos(np.pi * (x1 + 1) / (x1 - x2))
                    + 2 * x2**2
                    + 2 * np.pi**2 * x2
                    + np.pi**2
                ) / (4 * np.pi**2 * x1 - 4 * np.pi**2 * x2)
            elif abs(x1) <= 1 and abs(x2) <= 1:
                t1 = (4 * x1 + np.pi**2 * x1 - 4 * x2 + np.pi**2 * x2) / (4 * np.pi**2)
            elif abs(x1) <= 1 and x2 > 1:
                t1 = (
                    -2 * x1**2 * np.cos(np.pi * (x1 - 1) / (x1 - x2))
                    + 2 * x1**2
                    + np.pi**2 * x1**2
                    + 4 * x1 * x2 * np.cos(np.pi * (x1 - 1) / (x1 - x2))
                    - 4 * x1 * x2
                    - 2 * x2**2 * np.cos(np.pi * (x1 - 1) / (x1 - x2))
                    + 2 * x2**2
                    - 2 * np.pi**2 * x2
                    + np.pi**2
                ) / (4 * np.pi**2 * x1 - 4 * np.pi**2 * x2)
            elif x1 > 1 and x2 < -1:
                t1 = (
                    2
                    * x1**2
                    * np.sin(np.pi / (x1 - x2))
                    * np.sin(np.pi * x1 / (x1 - x2))
                    - 4
                    * x1
                    * x2
                    * np.sin(np.pi / (x1 - x2))
                    * np.sin(np.pi * x1 / (x1 - x2))
                    + np.pi**2 * x1
                    + 2
                    * x2**2
                    * np.sin(np.pi / (x1 - x2))
                    * np.sin(np.pi * x1 / (x1 - x2))
                    + np.pi**2 * x2
                ) / (2 * np.pi**2 * x1 - 2 * np.pi**2 * x2)
            elif x1 > 1 and abs(x2) <= 1:
                t1 = (
                    2 * x1**2 * np.cos(np.pi * (x1 - 1) / (x1 - x2))
                    + 2 * x1**2
                    - 4 * x1 * x2 * np.cos(np.pi * (x1 - 1) / (x1 - x2))
                    - 4 * x1 * x2
                    + 2 * np.pi**2 * x1
                    + 2 * x2**2 * np.cos(np.pi * (x1 - 1) / (x1 - x2))
                    - np.pi**2 * x2**2
                    + 2 * x2**2
                    - np.pi**2
                ) / (4 * np.pi**2 * x1 - 4 * np.pi**2 * x2)
            elif x1 > 1 and x2 > 1:
                t1 = 1 / 2

        y[n] = t0 + t1

        x2 = x1
        x1 = x0
    return y


if __name__ == "__main__":
    sampleRate = 48000

    frequencyHz = 1661
    durationSecond = 1
    gain = 100
    clipFunctions = [
        hardclipJ0,
        adaaCosineHardclipJ2,
        halfrectJ0,
        adaaCosineHalfRectJ2,
    ]

    for clipFn in clipFunctions:
        sig = generateSine(int(sampleRate * durationSecond), frequencyHz / sampleRate)
        # sig *= signal.windows.tukey(len(sig), 0.1 / durationSecond)
        sig = clipFn(sig, gain)
        sig = normalize(sig)

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
        ax[0].set_xlim([0, 3 * int(sampleRate / frequencyHz)])

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
