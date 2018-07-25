import matplotlib.pyplot as p
import numpy
import soundfile
import scipy.signal as signal

from pyfftw.interfaces.numpy_fft import rfft


def rcHighpass(samplerate, source, r, c):
    rc = r * c
    num, den, dt = signal.cont2discrete(([rc, 0], [rc, 1]), 1.0 / samplerate)
    return signal.lfilter(num[0], den, source)


def applyContinuousFilter(samplerate, source, system):
    num, den, dt = signal.cont2discrete(system, 1.0 / samplerate)
    return signal.lfilter(num[0], den, source)


def highMetalBandpass(samplerate, source):
    return applyContinuousFilter(samplerate, source, (
        [9.471e+10, 2.5625e+16, 0],
        [3437973.0, 1.816815e+11, 8.03e+15, 3.125e+20],
    ))


def lowMetalBandpass(samplerate, source):
    return applyContinuousFilter(samplerate, source, (
        [2.68345e+11, 3.5234375e+16, 0],
        [30108309.0, 5.227075e+11, 1.56671875e+16, 1.953125e+20],
    ))


def highMetal(samplerate, source):
    sig = highMetalBandpass(samplerate, source)
    sig = rcHighpass(samplerate, sig, 1e6, 47e-9)
    sig = rcHighpass(samplerate, sig, 5600, 470e-12)
    sig = rcHighpass(samplerate, sig, 100e3, 10e-10)
    return rcHighpass(samplerate, sig, 5600, 22e-10)


def lowMetal(samplerate, source):
    sig = lowMetalBandpass(samplerate, source)
    sig = rcHighpass(samplerate, sig, 1e6, 47e-9)
    sig = rcHighpass(samplerate, sig, 68000, 470e-12)
    return rcHighpass(samplerate, sig, 10000, 470e-12)


def toFrequencyResponse(data):
    dat = rfft(data, planner_effort="FFTW_ESTIMATE")
    return 20 * numpy.log10(numpy.abs(dat))


def getFreq(samplerate, response):
    freq_step = samplerate / 2 / (len(response) - 1)
    return numpy.arange(len(response)) * freq_step


samplerate = 44100

source = numpy.zeros(512)
source[0] = 1  # impulse

high_metal = highMetal(samplerate, source)
low_metal = lowMetal(samplerate, source)

res_high = toFrequencyResponse(high_metal)
res_low = toFrequencyResponse(low_metal)
x = getFreq(samplerate, res_high)

p.semilogx(x, res_high, label="High Metal")
p.semilogx(x, res_low, label="Low Metal")

p.xlabel("Frequency[Hz]")

p.ylabel("Gain[dB]")
p.ylim((-120, 40))

p.title("Frequency Response of Metal Filters")
p.grid(which="both")
p.legend()

p.show()
# p.savefig("metal_filter.png")
# p.close()
