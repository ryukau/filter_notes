import numpy as np
import scipy.signal as signal
import scipy.special as special
import matplotlib.pyplot as plt

def average(fir):
    denom = np.sum(fir)
    if denom == 0:
        return fir
    return fir / denom

def movingAverage(nTaps, integrate=0):
    length = nTaps // 2**integrate
    taps1 = average(np.ones(length))
    for _ in range(integrate):
        taps2 = average(np.ones(length + 1))  # フィルタのタップ数をそろえるために +1 。
        taps1 = average(signal.convolve(taps1, taps2))
        length *= 2
    return taps1

def thiranLowpassFixedGain(n, τ):
    a = np.empty(n + 1)
    for k in range(n + 1):
        i = np.arange(n + 1)
        a[k] = (-1)**k * special.binom(n, k) * np.prod((2 * τ + i) / (2 * τ + k + i))
    freq, resp = signal.freqz(1, a, worN=1)
    b = np.zeros_like(a)
    b[0] = 1 / resp.real
    return (b, a)

order = 4
delay = 512

fir_ma = movingAverage(delay, 1)

sos_bessel = signal.bessel(order, 2 / (np.pi * delay), norm="delay", output="sos")

b, a = thiranLowpassFixedGain(order, delay / 2)
sos_thiran = signal.tf2sos(b, a)

sig = np.hstack((np.zeros(delay), np.ones(delay * 2)))
out_trianglar = signal.lfilter(fir_ma, 1, sig)
out_bessel = signal.sosfilt(sos_bessel, sig)
out_thiran = signal.sosfilt(sos_thiran, sig)

plt.axvline(2 * delay, lw=1, color="black", alpha=0.5, ls="--")
plt.plot(sig, lw=1, color="black", label="Input")
plt.plot(out_trianglar, lw=1, alpha=0.75, color="red", label="Triangular")
plt.plot(out_bessel, lw=1, alpha=0.75, color="orange", label="Bessel")
plt.plot(out_thiran, lw=1, alpha=0.75, color="blue", label="Thiran")
plt.grid()
plt.legend()
plt.show()
