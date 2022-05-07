"""
Test DoubleAverageFilter and PeakHold.

DoubleAverageFilter output `smooth` and its reference do not match because of the
difference of float rounding mode.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile
from collections import deque

def peakHold(sig, holdTime, neutral=0):
    out = np.empty(len(sig))
    buffer = deque([0 for _ in range(holdTime)])
    hold = deque([])
    for i in range(len(sig)):
        x0 = sig[i]

        # Remove existing hold values which is less than new hold value.
        if len(hold) > 0:
            idx = len(hold) - 1
            while idx >= 0:
                if hold[idx] < x0:
                    hold.pop()
                else:
                    break
                idx -= 1

        hold.append(x0)  # Add new hold value to hold buffer.

        buffer.append(x0)
        delayOut = buffer.popleft()
        if len(hold) > 0 and delayOut == hold[0]:
            hold.popleft()

        out[i] = hold[0] if len(hold) > 0 else neutral

    return out

def getRefSmooth(refPeak, holdTime):
    fir = signal.get_window("bartlett", holdTime + 1)
    return signal.convolve(refPeak, fir / np.sum(fir))[:len(refPeak)]

def test(holdTime):
    sig, _ = soundfile.read("snd/input.wav")
    peak, _ = soundfile.read("output/input_peak.wav")
    smooth, _ = soundfile.read("output/input_smooth.wav")

    delayed = np.hstack((np.zeros(holdTime), sig[:-holdTime]))
    refPeak = peakHold(sig, holdTime)
    refSmooth = getRefSmooth(refPeak, holdTime)

    eps = np.finfo(np.float32).eps
    np.testing.assert_allclose(smooth, refSmooth, rtol=1e-3)
    np.testing.assert_allclose(peak, refPeak, rtol=eps)

    plt.figure()
    plt.plot(delayed, lw=1, color="orange", alpha=0.5, label="Input")
    plt.plot(smooth, lw=1, color="red", alpha=0.5, label="C++")
    plt.plot(refSmooth, lw=1, color="blue", alpha=0.5, label="target")
    plt.grid()
    plt.legend()
    # plt.savefig("img/DoubleAverageFilterTest.svg")

    plt.figure()
    plt.plot(sig, lw=1, color="orange", alpha=0.5, label="sig")
    plt.plot(peak, lw=1, color="red", alpha=0.5, label="peak")
    plt.plot(refPeak, lw=1, color="blue", alpha=0.5, label="target")
    plt.grid()
    plt.legend()
    plt.show()

if __name__ == "__main__":
    test(32)
