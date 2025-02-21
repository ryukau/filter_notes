import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.signal as signal
import soundfile
from collections import deque


def adaa1Bypass(x):
    eps = np.finfo(np.float64).eps
    y = np.zeros_like(x)
    x1 = 0
    s1 = 0
    for n in range(len(x)):
        x0 = x[n]
        s0 = x0 * x0 / 2
        if (x1 == 0 and s1 == 0) or np.abs(x0 - x1) < eps:
            y[n] = (x0 + x1) / 2
        else:
            y[n] = (s0 - s1) / (x0 - x1)
        s1 = s0
        x1 = x0
    return y


def adaa2Bypass(x):
    def f1(x):
        return x * x / 2

    def f2(x):
        return x * x * x / 6

    eps = np.finfo(np.float64).eps
    y = np.zeros_like(x)

    x1 = 0
    x2 = 0
    s1 = 0
    for n in range(len(x)):
        x0 = x[n]

        f2_x1 = f2(x1)
        s0 = (
            f1((x0 + x1) / 2) if np.abs(x0 - x1) < eps else (f2(x0) - f2_x1) / (x0 - x1)
        )

        if x1 == 0 and x2 == 0:
            y[n] = (x0 + 2 * x1 + x2) / 4
        elif np.abs(x0 - x2) < eps:
            x_bar = (x0 + x2) / 2
            delta = x_bar - x1
            if np.abs(delta) < eps:
                y[n] = (x_bar + x1) / 2
            else:
                y[n] = (2 / delta) * (f1(x_bar) + (f2_x1 - f2(x_bar)) / delta)
        else:
            y[n] = 2 * (s0 - s1) / (x0 - x2)
        s1 = s0
        x2 = x1
        x1 = x0
    return y


def adaa2CosineWindow(x):
    eps = np.finfo(np.float64).eps
    x = x.copy()
    y = np.zeros_like(x)

    x1 = 0
    x2 = 0
    for n, _ in enumerate(x):
        x0 = x[n]

        d0 = x0 - x1
        t0 = 0
        if np.abs(d0) < eps:
            mid = 0.5 + 2 / (np.pi * np.pi)
            t0 = 0.5 * (x2 + mid * (x1 - x2))
        else:
            t0 = (-x0 + x1 + np.pi**2 * (x0 + x1) / 4) / np.pi**2

        d1 = x1 - x2
        t1 = 0
        if np.abs(d1) < eps:
            mid = 0.5 - 2 / (np.pi * np.pi)
            t1 = 0.5 * (x1 + mid * (x2 - x1))
        else:
            t1 = (x1 - x2 + np.pi**2 * (x1 + x2) / 4) / np.pi**2

        y[n] = t0 + t1

        x2 = x1
        x1 = x0
    return y


def peakHold(sig, holdTime):
    """
    sig: 入力信号。
    holdTime: ホールド時間。単位はサンプル数。
    """
    out = np.empty(len(sig))
    buffer = deque([0 for _ in range(holdTime)])
    hold = deque([])
    for i in range(len(sig)):
        x0 = sig[i]

        if len(hold) > 0:
            idx = len(hold) - 1
            while idx >= 0:
                if hold[idx] < x0:
                    hold.pop()
                else:
                    break
                idx -= 1

        hold.append(x0)

        buffer.append(x0)
        delayOut = buffer.popleft()

        if len(hold) > 0 and delayOut == hold[0]:
            hold.popleft()

        out[i] = hold[0] if len(hold) > 0 else 0

    return out


def applyCharacteristicCurve(amp, threshold):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(amp > threshold, threshold / amp, 1)


def getTriangleWindow(length):
    win = signal.get_window("bartlett", length)[::-1]
    win /= np.sum(win)
    return win


def applyLimiter(sig, absed, holdTime):
    hold = peakHold(absed, holdTime)
    gain = applyCharacteristicCurve(hold, 1.0)
    fir = getTriangleWindow(holdTime)
    smoothed = np.convolve(gain, fir)[holdTime - 1 :]
    return sig * smoothed


def generateSignal(upRate, length, seed=None):
    rng = np.random.default_rng(seed)
    sig = 2.0 * rng.binomial(1, 0.5, length) - 1

    # frequencyNormalized = 16661 / 48000
    # phase = 0
    # sig = np.empty(length)
    # for i in range(len(sig)):
    #     sig[i] = np.sin(2 * np.pi * phase)
    #     phase += frequencyNormalized
    #     phase -= np.floor(phase)

    up = signal.resample(sig, upRate * len(sig))
    downAbs = np.max(np.abs(up).reshape((len(sig), upRate)), axis=1)
    return (sig, up, downAbs)


def testGenerateSignal():
    length = 16
    seed = 6871867313
    sig, up, downAbs = generateSignal(16, length, seed)

    plt.plot(np.abs(up))
    plt.plot(np.repeat(downAbs, 16))
    plt.grid()
    plt.show()


if __name__ == "__main__":
    length = 65536
    holdTime = 128
    upRate = 16
    seed = None

    sig, up, downAbs = generateSignal(upRate, length, seed)

    src = sig.copy()
    # sig = adaa1Bypass(sig)
    # sig = adaa2Bypass(sig)
    sig = adaa2CosineWindow(sig)
    limited = applyLimiter(sig, np.abs(sig), holdTime)
    # limited = applyLimiter(sig, downAbs, holdTime) # Reference quality.

    soundfile.write("src.wav", src, 48000, "FLOAT")
    soundfile.write("out.wav", sig, 48000, "FLOAT")

    upLimited = signal.resample(limited, upRate * len(limited))
    peak = np.max(np.abs(upLimited))
    print(f"Peak [dB]: {20 * math.log10(peak)}")

    peakIndex = np.where(np.abs(upLimited) > 1)[0]
    plt.plot(up, lw=1, color="red", alpha=0.5, label="16x Source")
    plt.plot(upLimited, lw=1, color="black", label="16x Limited")
    plt.legend(loc="lower right")
    plt.grid(color="#f0f0f0")
    plt.tight_layout()
    plt.show()
