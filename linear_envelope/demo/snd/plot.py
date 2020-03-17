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

adsr, samplerate = soundfile.read("ADSR.wav", always_2d=True)
plot("ADSR", adsr.T[0], samplerate)

relAtk, samplerate = soundfile.read("ReleaseWhileAttack.wav", always_2d=True)
plot("ReleaseWhileAttack", relAtk.T[0], samplerate)

relDec, samplerate = soundfile.read("ReleaseWhileDecay.wav", always_2d=True)
plot("ReleaseWhileDecay", relDec.T[0], samplerate)

trigRel, samplerate = soundfile.read("TriggerWhileRelease.wav", always_2d=True)
plot("TriggerWhileRelease", trigRel.T[0], samplerate)

pCtrl, samplerate = soundfile.read("PController.wav", always_2d=True)
plot("PController", pCtrl.T[0], samplerate)
