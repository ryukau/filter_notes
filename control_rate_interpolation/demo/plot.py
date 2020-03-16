import matplotlib.pyplot as pyplot
import numpy
import soundfile

def plot(name: str, data, samplerate: float):
    xTick = numpy.arange(len(data)) / samplerate
    pyplot.figure(figsize=(8, 2))
    pyplot.title(name)
    pyplot.plot(xTick, data)
    pyplot.ylim([-0.1, 1.1])
    pyplot.grid()
    pyplot.show()

def compare(name1: str, data1, name2: str, data2, samplerate: float):
    xTick = numpy.arange(len(data1)) / samplerate
    pyplot.figure(figsize=(8, 2))
    pyplot.plot(xTick, data1, color="red", lw=1, alpha=0.75, label=name1)
    pyplot.plot(xTick, data2, color="blue", lw=1, alpha=0.75, label=name2)
    # pyplot.ylim([-0.1, 1.1])
    pyplot.legend()
    pyplot.grid()
    pyplot.show()

def writeSin(name, samplerate, envelope):
    frequency = 100
    time = len(envelope) / samplerate
    phase = numpy.linspace(0, 2 * numpy.pi * frequency * time, len(envelope))
    sin = numpy.sin(phase)
    wav = sin * envelope
    soundfile.write(f"snd/sin_{name}.wav", wav, samplerate, subtype="FLOAT")

def processExample(name):
    sig, fs = soundfile.read(f"snd/{name}.wav", always_2d=True)
    sig = sig.T[0]
    plot(name, sig, fs)
    writeSin(name, fs, sig)

processExample("naive")
processExample("linterp")
processExample("smoother")

exit()

# bench_linterp.cpp
linterpRamp, samplerate = soundfile.read("snd/linterpRamp.wav", always_2d=True)
linterpIndex, samplerate = soundfile.read("snd/linterpIndex.wav", always_2d=True)
compare("linterpRamp", linterpRamp.T[0], "linterpIndex", linterpIndex.T[0], samplerate)

# bench_smoother.cpp
smootherRamp, samplerate = soundfile.read("snd/smootherRamp.wav", always_2d=True)
smootherIndex, samplerate = soundfile.read("snd/smootherIndex.wav", always_2d=True)
compare(
    "smootherRamp",
    smootherRamp.T[0],
    "smootherIndex",
    smootherIndex.T[0],
    samplerate,
)

smootherRampCommon, samplerate = soundfile.read(
    "snd/smootherRampCommon.wav", always_2d=True)
smootherIndexCommon, samplerate = soundfile.read(
    "snd/smootherIndexCommon.wav", always_2d=True)
compare(
    "smootherRampCommon",
    smootherRampCommon.T[0],
    "smootherIndexCommon",
    smootherIndexCommon.T[0],
    samplerate,
)
