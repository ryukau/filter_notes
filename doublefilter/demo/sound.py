import dent
import soundfile
import numpy as np
import scipy.signal as signal

import matplotlib.pyplot as plt

class Model:
    def __init__(self):
        self.acc1 = 0
        self.vel1 = 0
        self.pos1 = 0

        self.acc2 = 0
        self.vel2 = 0
        self.pos2 = 0

        self.x1 = 0

    def process(self, k1, k2, altGain, isHighpass, x0):
        v2Gain = np.sqrt(k1) if altGain else 1

        self.acc2 = k2 * (self.vel1 - self.vel2)
        self.vel2 += self.acc2 + x0 - self.x1
        self.pos2 += self.vel2 * v2Gain

        self.acc1 = -k1 * self.pos1 - self.acc2
        self.vel1 += self.acc1
        self.pos1 += self.vel1

        self.x1 = x0

        if isHighpass:
            self.pos1 *= 0.999
            return self.pos1
        self.pos2 *= 0.999
        return self.pos2

def cutoffCurve(samplerate, cutoffHz):
    x = cutoffHz / samplerate
    return 6.5451144600705975 * x + 20.46391326872472 * x * x

def normalize(sig):
    m = np.max(np.abs(sig))
    if m == 0:
        return sig
    return sig / m

samplerate = 48000
duration = 1
frequency = 45
nSample = int(duration * samplerate)

time = np.linspace(0, frequency * duration, nSample)
phase = 2 * np.pi * time
sig = 0.01 * signal.sawtooth(phase)

def render(sig, k2Array, prefix, resonance, altGain, isHighpass):
    _, k1Array = dent.resonanceCurve(resonance, k2Array,
                                     dent.dentK1 if altGain else dent.shelveK1)
    model = Model()
    output = [
        model.process(k1, k2, altGain, isHighpass, x0)
        for k1, k2, x0 in zip(k1Array, k2Array, sig)
    ]
    soundfile.write(
        f"{prefix}_res{resonance:.5f}_altGain{altGain}_isHighpass{isHighpass}.wav",
        normalize(output),
        samplerate,
        subtype="FLOAT")

cutoff = 5000 * np.geomspace(1e-5, 1, nSample)
k2Array = cutoffCurve(samplerate, cutoff)

render(sig, k2Array, "exp", 1, False, False)
render(sig, k2Array, "exp", 0.1, False, False)
render(sig, k2Array, "exp", 0.001, False, False)

render(sig, k2Array, "exp", 1, True, False)
render(sig, k2Array, "exp", 0.1, True, False)
render(sig, k2Array, "exp", 0.001, True, False)

render(sig, k2Array, "exp", 1, False, True)
render(sig, k2Array, "exp", 0.1, False, True)
render(sig, k2Array, "exp", 0.001, False, True)

render(sig, k2Array, "exp", 1, True, True)
render(sig, k2Array, "exp", 0.1, True, True)
render(sig, k2Array, "exp", 0.001, True, True)

cutoff = 5000 * np.linspace(1e-5, 1, nSample)
k2Array = cutoffCurve(samplerate, cutoff)

render(sig, k2Array, "lin", 1, False, False)
render(sig, k2Array, "lin", 0.1, False, False)
render(sig, k2Array, "lin", 0.001, False, False)

render(sig, k2Array, "lin", 1, True, False)
render(sig, k2Array, "lin", 0.1, True, False)
render(sig, k2Array, "lin", 0.001, True, False)

render(sig, k2Array, "lin", 1, False, True)
render(sig, k2Array, "lin", 0.1, False, True)
render(sig, k2Array, "lin", 0.001, False, True)

render(sig, k2Array, "lin", 1, True, True)
render(sig, k2Array, "lin", 0.1, True, True)
render(sig, k2Array, "lin", 0.001, True, True)
