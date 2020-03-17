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

def process(name):
    data, samplerate = soundfile.read(f"snd/{name}.wav", always_2d=True)
    data = data.T[0]
    plot(name, data, samplerate)

    frequency = 100
    phase = numpy.linspace(0, 2 * numpy.pi * frequency * len(data) / samplerate,
                           len(data))
    tone = numpy.sin(phase)
    soundfile.write(f"snd/tone_{name}.wav", data * tone, samplerate, subtype="FLOAT")

process("ADSR")
process("ReleaseWhileAttack")
process("ReleaseWhileDecay")
process("TriggerWhileRelease")
process("ChangeSustain")
