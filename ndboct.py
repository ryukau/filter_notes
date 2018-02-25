import math
import matplotlib.pyplot as p
import numpy as np
from scipy import signal
from scipy.fftpack import *


def showResponse(sample_rate, sig):
    res = fftshift(fft(sig, 2**16) / (len(sig) / 2.0))
    gain = np.array_split(
        20 * np.log10(abs(res / np.max(abs(res)))),
        2,
    )[0][::-1]
    freq = getFreq(sample_rate, gain)
    p.plot(freq, gain)
    p.grid(which='both')
    p.show()


def makeFrequencySpecification(sample_rate, length, low, high, slope):
    # slope は負の値。単位は[dB/oct]。
    # low, high は周波数[Hz]。
    #
    # dB/oct を算出する式の解説。
    # https://femci.gsfc.nasa.gov/random/randomequations.html
    #
    # 6.020599913279624 = 20 * np.log10(2)
    #
    oct = (high / low)**(slope / 6.020599913279624)
    spec = [0] * length
    freq_step = sample_rate / 2 / (length - 1)
    for i in range(0, len(spec)):
        freq = i * freq_step
        if freq <= low:
            spec[i] = 1
        elif freq >= high:
            spec[i] = oct
        else:
            spec[i] = (freq / low)**(slope / 6.020599913279624)
    return spec


def getFreq(sample_rate, spec):
    freq_step = sample_rate / 2 / (len(spec) - 1)
    return [i * freq_step for i in range(0, len(spec))]


def toDecibel(spec):
    return [20 * math.log10(v) for v in spec]


if __name__ == "__main__":
    sample_rate = 44100
    spec = makeFrequencySpecification(sample_rate, 1024, 100, 10000, -10)

    # 設計した特性。
    spec_db = toDecibel(spec)
    freq = getFreq(sample_rate, spec)

    p.title('Specification')
    p.plot(freq, spec_db)
    p.xscale('log')
    p.grid(which='both')
    p.show()

    # フィルタ係数。
    filt = fftshift(ifft(spec + spec[::-1]))

    p.title('Filter Coefficients')
    p.plot(filt)
    p.show()

    p.title('Frequency Response of the filter (Log)')
    p.xscale('log')
    showResponse(sample_rate, filt)

    # x軸が線形の場合をプロット。
    p.title('Frequency Response of the filter (Linear)')
    showResponse(sample_rate, filt)

    # フィルタを乱数ノイズに適用。
    noise = np.random.uniform(-1.0, 1.0, sample_rate)

    p.title('Before')
    p.xscale('log')
    showResponse(sample_rate, noise)

    filtered = signal.convolve(noise, filt, mode='same')

    p.title('After')
    p.xscale('log')
    showResponse(sample_rate, filtered)
