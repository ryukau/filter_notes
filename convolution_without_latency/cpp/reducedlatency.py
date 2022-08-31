import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import soundfile
import json

def firlp2hp(lp):
    """
    `index` search is heuristic. Check if something is wrong.
    """
    hp = -lp
    index = np.argmax(np.abs(hp))
    hp[index] += 1
    return hp

def generateSin(frequency, duration, sampleRate=48000):
    length = int(sampleRate * duration)
    phase = np.linspace(0, 2 * np.pi * frequency * duration, length)
    return (np.sin(phase), sampleRate)

def overlapAddNaive(sig, coefficient, fs):
    fftsize = 2 * len(coefficient)
    half = fftsize // 2

    fir = np.fft.rfft(coefficient, fftsize)

    source = np.hstack([sig, np.zeros(fftsize - len(sig) % fftsize)])
    output = np.zeros_like(source)
    buf = np.zeros(fftsize)

    index: int = 0
    while index < len(source) - half:
        buf[:half] = source[index:index + half]

        spc = np.fft.rfft(buf)
        spc *= fir
        flt = np.fft.irfft(spc)

        output[index:index + fftsize] += flt

        index += half
    return output

def naiveConvolve(sigA, sigB):
    if len(sigA) < len(sigB):
        sigA, sigB = (sigB, sigA)

    padded = np.hstack([np.zeros(len(sigB) - 2), sigA, [0]])

    output = np.zeros(len(padded))
    for idx in range(len(output)):
        for jdx in range(len(sigB)):
            if idx + jdx >= len(padded):
                break
            output[idx] += sigB[jdx] * padded[idx + jdx]
    return output

def overlapSave(sig, coefficient, fs):
    fftsize = 2 * len(coefficient)
    half = fftsize // 2

    fir = np.fft.rfft(coefficient, fftsize)

    source = np.hstack([sig, np.zeros(fftsize - len(sig) % fftsize)])
    output = np.zeros_like(source)
    buf = np.zeros(fftsize)

    index: int = 0
    while index < len(source) - half:
        buf[:half] = buf[half:]
        buf[half:] = source[index:index + half]

        spc = np.fft.rfft(buf)
        spc *= fir
        flt = np.fft.irfft(spc)

        output[index:index + half] += flt[half:]

        index += half
    return output

def minimumCost(sig, coefficient, fs):
    fftsize = 2 * len(coefficient)
    half = fftsize // 2

    segment = len(coefficient) // 2
    h0 = coefficient[:segment]
    h1 = coefficient[segment:]

    fir0 = np.fft.rfft(h0, 2 * len(h0))
    fir1 = np.fft.rfft(h1, 2 * len(h1))

    source = np.hstack([sig, np.zeros(fftsize - len(sig) % fftsize)])
    output = np.zeros_like(source)

    buf0 = np.zeros(2 * len(h0))
    buf1 = np.zeros(2 * len(h1))

    half0 = len(buf0) // 2
    half1 = len(buf1) // 2

    index: int = 0
    while index < len(source) - half0:
        buf0[:half0] = buf0[half0:]
        buf0[half0:] = source[index:index + half0]

        spc0 = np.fft.rfft(buf0)
        spc0 *= fir0
        flt0 = np.fft.irfft(spc0)

        output[index:index + half0] += flt0[half0:]

        index += half0

    index = 0
    while index < len(source) - 2 * half1:
        buf1[:half1] = buf1[half1:]
        buf1[half1:] = source[index:index + half1]

        spc1 = np.fft.rfft(buf1)
        spc1 *= fir1
        flt1 = np.fft.irfft(spc1)

        output[index + half1:index + 2 * half1] += flt1[half1:]

        index += half1

    return output

coefficient = signal.firwin(2047, 1000, window="nuttall", fs=48000)
coefficient = np.hstack([firlp2hp(coefficient), [0]])

sigA, fs = generateSin(1000, 0.1)
sigB, fs = generateSin(100, 0.1)
sig = sigA + sigB

# fs = 48000
# sig = np.hstack([np.zeros(len(coefficient) - 1), np.ones(int(1.0 * fs))])

add = overlapAddNaive(sig, coefficient, fs)
# save = overlapSave(sig, coefficient, fs)
minCost = minimumCost(sig, coefficient, fs)
scipy_convolve = signal.convolve(sig, coefficient)
# naive = naiveConvolve(sig, coefficient)

# plt.plot(sig, alpha=0.5, label="input")
# plt.plot(add, alpha=0.5, label="add")
# plt.plot(save, alpha=0.5, label="save")
plt.plot(minCost, alpha=0.5, label="minCost")
plt.plot(scipy_convolve, alpha=0.5, label="SciPy")
# plt.plot(naive, alpha=0.5, label="naive")
plt.grid()
plt.legend()
plt.show()
