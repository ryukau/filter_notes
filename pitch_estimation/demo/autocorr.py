import matplotlib.pyplot as pyplot
import numpy
import scipy.signal
from pyfftw.interfaces.numpy_fft import fft, ifft

def generate_sin(duration, samplerate, frequency):
    length = int(duration * samplerate)
    phase = numpy.linspace(0, 2 * numpy.pi * frequency * duration, length)
    return numpy.sin(phase)

def autocorrelation_naive(sig):
    corr = numpy.zeros(len(sig))
    for r in range(len(sig)):
        for k in range(len(sig) - r):
            corr[r] += sig[k] * sig[k + r]
    return corr

def autocorr_fft(sig):
    len_sig = len(sig)
    sig = numpy.pad(sig, (0, len_sig), "constant")
    spec = fft(sig)
    return ifft(spec * spec.conj()).real[:len_sig]

def autocorr_numpy(sig):
    corr = numpy.correlate(sig, sig, mode="full")
    return corr[len(sig) - 1::-1]

sin_wave = generate_sin(0.1, 44100, 100)

corr_naive = autocorrelation_naive(sin_wave)
corr_fft = autocorr_fft(sin_wave)
corr_numpy = autocorr_numpy(sin_wave)

numpy.testing.assert_allclose(corr_naive, corr_numpy, atol=1e-5)
numpy.testing.assert_allclose(corr_fft, corr_numpy, atol=1e-5)

cmap = pyplot.get_cmap("plasma")
pyplot.plot(corr_naive, label="naive", lw=10, color=cmap(0.0))
pyplot.plot(corr_fft, label="fft", lw=6, color=cmap(0.5))
pyplot.plot(corr_numpy, label="numpy", lw=2, color=cmap(1.0))
pyplot.legend()
pyplot.show()
