import numpy
import math
import matplotlib.pyplot as pyplot
import soundfile
from scipy.signal import windows

def to_decibel(data):
    data_abs = numpy.abs(data)
    return 20 * numpy.log10(data_abs / numpy.max(data_abs))

def signal_to_spectrum(sig):
    return to_decibel(numpy.abs(numpy.fft.rfft(sig)))

def blamp_residual4(d):
    return (
        d**5 / 120,
        -d**5 / 40 + d**4 / 24 + d**3 / 12 + d**2 / 12 + d / 24 + 1 / 120,
        d**5 / 40 - d**4 / 12 + d**2 / 3 - d / 2 + 7 / 30,
        -d**5 / 120 + d**4 / 24 - d**3 / 12 + d**2 / 12 - d / 24 + 1 / 120,
    )

# yapf: disable
def blamp_residual6(t):
    return numpy.array((
        t**7/5040,
        -t**7/1008+t**6/720+t**5/240+t**4/144+t**3/144+t**2/240+t/720+1/5040,
        t**7/504-t**6/180-t**5/120+t**4/72+(5*t**3)/72+(13*t**2)/120+(29*t)/360+61/2520,
        -t**7/504+t**6/120-t**4/24+(11*t**2)/40-t/2+239/840,
        t**7/1008-t**6/180+t**5/120+t**4/72-(5*t**3)/72+(13*t**2)/120-(29*t)/360+61/2520,
        -t**7/5040+t**6/720-t**5/240+t**4/144-t**3/144+t**2/240-t/720+1/5040,
    ))

def blamp_residual8(t):
    return numpy.array((
        t**9/362880,
        -t**9/51840+t**8/40320+t**7/10080+t**6/4320+t**5/2880+t**4/2880+t**3/4320+t**2/10080+t/40320+1/362880,
        t**9/17280-t**8/6720-t**7/2520+t**5/360+t**4/120+(7*t**3)/540+t**2/84+(31*t)/5040+1/720,
        -t**9/10368+t**8/2688+t**7/2016-t**6/480-(19*t**5)/2880+t**4/192+(49*t**3)/864+(397*t**2)/3360+(4541*t)/40320+347/8064,
        t**9/10368-t**8/2016+t**6/270-t**4/36+(151*t**2)/630-t/2+1487/4536,
        -t**9/17280+t**8/2688-t**7/2016-t**6/480+(19*t**5)/2880+t**4/192-(49*t**3)/864+(397*t**2)/3360-(4541*t)/40320+347/8064,
        t**9/51840-t**8/6720+t**7/2520-t**5/360+t**4/120-(7*t**3)/540+t**2/84-(31*t)/5040+1/720,
        -t**9/362880+t**8/40320-t**7/10080+t**6/4320-t**5/2880+t**4/2880-t**3/4320+t**2/10080-t/40320+1/362880,
    ))
# yapf: enable

class TriangleNaive:
    def __init__(self, samplerate, frequency):
        self.phase = 0
        self.tick = frequency / samplerate

    def process(self):
        output = 4 * numpy.abs(self.phase - 0.5) - 1
        self.phase += self.tick
        if self.phase > 1.0:
            self.phase -= 1.0
        return output

class TriangleBLAMP:
    def __init__(self, points, samplerate, frequency, residual_func):
        self.residual_func = residual_func

        self.s = numpy.zeros(points)
        self.isEdge = numpy.zeros(points)

        self.edgePoint = points // 2 - 1

        self.phase = 0
        self.tick = frequency / samplerate
        self.mu = 4 * frequency / samplerate
        self.d = 0

    def process(self):
        self.s[-1] = 4 * numpy.abs(self.phase - 0.5) - 1

        if self.isEdge[self.edgePoint] == 1:
            self.s += self.mu * (
                self.residual_func(self.d) + numpy.flip(self.residual_func(1 - self.d)))
        elif self.isEdge[self.edgePoint] == -1:
            self.s -= self.mu * (
                self.residual_func(self.d) + numpy.flip(self.residual_func(1 - self.d)))

        self.phase += self.tick
        if self.phase >= 0.5 and self.phase < 0.5 + self.tick:
            self.isEdge[-1] = 1
            self.d = (self.phase - 0.5) / self.tick
        elif self.phase > 1.0:
            self.phase -= 1.0
            self.isEdge[-1] = -1
            self.d = self.phase / self.tick
        else:
            self.isEdge[-1] = 0

        self.s = numpy.roll(self.s, -1)
        self.isEdge = numpy.roll(self.isEdge, -1)
        return self.s[0]

def render_blamp(points, samplerate, frequency, residual_func):
    triangle = TriangleBLAMP(points, samplerate, frequency, residual_func)
    data = numpy.empty(n_sample)
    for _ in range(points - 2):
        triangle.process()
    for i in range(len(data)):
        data[i] = triangle.process()
    return data

def render_naive(samplerate, frequency):
    triangle = TriangleNaive(samplerate, frequency)
    data = numpy.empty(n_sample)
    for i in range(len(data)):
        data[i] = triangle.process()
    return data

def plot_waveform(data):
    """
    data = [{"name": "Tri0", "wave": tri0_wave, "spectrum": tri0_spec}, ...]
    """
    fig, ax = pyplot.subplots(len(data))
    fig.set_size_inches(6.4, 8.0)
    fig.tight_layout(pad=2)

    for i, item in enumerate(data):
        ax[i].set_title(item["name"])
        ax[i].plot(item["wave"][0:256], lw=1, color="black")
        ax[i].set_ylim((-1.1, 1.1))
        ax[i].grid()
    pyplot.show()

def plot_spectrum(data):
    fig, ax = pyplot.subplots(len(data))
    fig.set_size_inches(6.4, 8.0)
    fig.tight_layout(pad=2)

    for i, item in enumerate(data):
        ax[i].set_title(item["name"])
        ax[i].plot(item["spectrum"], lw=1, color="black")
        ax[i].set_ylim((-100, 0))
        ax[i].grid()
    pyplot.show()

def write_wav(samplerate, data):
    for item in data:
        soundfile.write(
            f"snd/{item['name']}.wav", item["wave"], samplerate, subtype="FLOAT")

samplerate = 44100
frequency = 5500
duration = 1.0
n_sample = int(duration * samplerate)

tri_naive = render_naive(samplerate, frequency)
tri_blamp4 = render_blamp(4, samplerate, frequency, blamp_residual4)
tri_blamp6 = render_blamp(6, samplerate, frequency, blamp_residual6)
tri_blamp8 = render_blamp(8, samplerate, frequency, blamp_residual8)

data = [
    {
        "name": "tri_naive",
        "wave": tri_naive,
        "spectrum": signal_to_spectrum(tri_naive),
    },
    {
        "name": "tri_blamp4",
        "wave": tri_blamp4,
        "spectrum": signal_to_spectrum(tri_blamp4),
    },
    {
        "name": "tri_blamp6",
        "wave": tri_blamp6,
        "spectrum": signal_to_spectrum(tri_blamp6),
    },
    {
        "name": "tri_blamp8",
        "wave": tri_blamp8,
        "spectrum": signal_to_spectrum(tri_blamp8),
    },
]

plot_waveform(data)
plot_spectrum(data)
write_wav(samplerate, data)
