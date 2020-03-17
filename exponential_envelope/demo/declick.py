import numpy
import matplotlib.pyplot as pyplot
import soundfile

samplerate = 48000
duration = 0.1
declick = 0.001

curve = numpy.hstack((
    numpy.geomspace(1, 1e-5, int(duration * samplerate)),
    numpy.zeros(samplerate),
))

phase = numpy.linspace(-numpy.pi, 0, int(declick * samplerate))
subEnv = numpy.hstack((
    0.5 + 0.5 * numpy.cos(phase),
    numpy.ones(len(curve) - len(phase)),
))
declicked = curve * subEnv

pyplot.figure(figsize=(6, 3))
pyplot.title("Declick")
pyplot.plot(declicked, alpha=0.5, color="blue", label="Declicked")
pyplot.plot(curve, alpha=0.5, color="red", label="Raw")
pyplot.xlabel("Samples")
pyplot.ylabel("Amplitude")
pyplot.xlim([-10, 200])
pyplot.grid()
pyplot.legend()
pyplot.tight_layout()
pyplot.show()

frequency = 1000
phase = numpy.linspace(0, 2 * numpy.pi * frequency * len(curve) / samplerate, len(curve))
tone = numpy.sin(phase + numpy.pi / 2)
soundfile.write(f"snd/declick_on.wav", declicked * tone, samplerate, subtype="FLOAT")
soundfile.write(f"snd/declick_off.wav", curve * tone, samplerate, subtype="FLOAT")
