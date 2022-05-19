import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.linalg as linalg
import scipy.signal as signal
import json
from collections import deque

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
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(amp > threshold, threshold / amp, 1)

def getTriangleWindow(length):
    win = signal.get_window("bartlett", length)[::-1]
    win /= np.sum(win)
    return win

def applyLimiter(sig, absed, holdTime):
    hold = peakHold(absed, holdTime)
    gain = applyCharacteristicCurve(hold, 1.0)
    fir = getTriangleWindow(holdTime)
    smoothed = np.convolve(gain, fir)[holdTime - 1:]
    return sig * smoothed

def generateSignal(upRate, seed=None):
    rng = np.random.default_rng(seed)
    sig = 2.0 * rng.binomial(1, 0.5, length) - 1
    up = signal.resample(sig, upRate * len(sig))
    downAbs = np.max(np.abs(up).reshape((len(sig), upRate)), axis=1)
    return (sig, up, downAbs)

def testGenerateSignal():
    length = 16
    seed = 6871867313
    sig, up, downAbs = generateSignal(seed)

    plt.plot(np.abs(up))
    plt.plot(np.repeat(downAbs, 16))
    plt.grid()
    plt.show()

if __name__ == "__main__":
    length = 2048
    holdTime = 128
    upRate = 16
    seed = 6415456
    fix = 10**(-0.2 / 20)  # Safe guard gain. Not used for demonstration.

    sig, up, downAbs = generateSignal(upRate, seed)

    # This implementation skips pre-lowpass which applies to `sig` to reduce overshoots.

    limited = applyLimiter(sig, downAbs, holdTime)

    upLimited = signal.resample(limited, upRate * len(limited))
    peak = np.max(np.abs(upLimited))
    print(f"Peak [dB]: {20 * math.log10(peak)}")
    # upLimited *= fix  # Safe guard gain is disabled for demonstration.

    peakIndex = np.where(np.abs(upLimited) > 1)[0]
    plt.plot(up, lw=1, color="red", alpha=0.5, label="16x Source")
    plt.plot(upLimited, lw=1, color="black", label="16x Limited")
    plt.scatter(peakIndex, upLimited[peakIndex], marker="s", zorder=10, label="Overshoot")
    plt.legend(loc="center right")
    plt.grid(color="#f0f0f0")
    plt.tight_layout()
    plt.show()
