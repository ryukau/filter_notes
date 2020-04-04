import numpy as np
import scipy.signal as signal
import math
import soundfile

class LP3:
    def __init__(self):
        self.reset()

    def reset(self):
        self.acc = 0
        self.vel = 0
        self.pos = 0
        self.x1 = 0

    def setK(self, c, resonance, uniformPeak):
        if not uniformPeak:
            return max(0, min(resonance, 0.99999))

        kExp = math.exp(-5.6852537097945195 * resonance)
        kMin = 1 - kExp
        kMax = 0.9999771732485103 - 0.01 * (kExp - 0.0033956716251850594)
        return kMax - (kMax - kMin) * math.acos(1 - c) / (np.pi / 2)

    def process(
        self,
        samplerate,
        lowpassHz,
        highpassHz,
        resonance,
        uniformPeak,
        uniformGain,
        x0,
    ):
        x = lowpassHz / samplerate
        c = (56.85341479156533 * x * x * x * x * x * x -
             60.92051508862034 * x * x * x * x * x - 1.6515635438744682 * x * x * x * x +
             31.558896956675998 * x * x * x - 20.61402812645397 * x * x +
             6.320753515093109 * x)

        k = self.setK(c, resonance, uniformPeak)

        y = highpassHz / samplerate
        alpha = (0.5638865655409118 +
                 0.43611343445908823 * math.exp(-6.501239408777854 * y))

        self.acc = c * self.vel + k * self.acc
        self.vel -= self.acc + x0 - self.x1
        self.pos = alpha * (self.pos - (c / (1 - k) if uniformGain else c) * self.vel)

        self.x1 = x0
        return self.pos

def normalize(sig):
    m = np.max(np.abs(sig))
    if m == 0:
        return sig
    return sig / m

def render(prefix, lp3, sig, cutoff, resonance, uniformPeak):
    lp3.reset()
    filtered = np.array([
        lp3.process(samplerate, low, 20, resonance, uniformPeak, False, x)
        for low, x in zip(cutoff, sig)
    ])
    soundfile.write(
        f"{prefix}_uniformPeak{uniformPeak}_resonance{resonance:.1f}.wav",
        normalize(filtered),
        samplerate,
        subtype="FLOAT")

samplerate = 48000
duration = 0.2
frequency = 45
section = 32

maxCutoff = np.geomspace(10, 22000, section)
cutoff = np.hstack(
    [peak * np.geomspace(1, 1e-5, int(duration * samplerate)) for peak in maxCutoff])

time = np.linspace(0, frequency * duration * section, len(cutoff))
phase = 2 * np.pi * time
sig = 0.01 * signal.sawtooth(phase)

lp3 = LP3()

render("exp", lp3, sig, cutoff, 0, False)
render("exp", lp3, sig, cutoff, 0.5, False)
render("exp", lp3, sig, cutoff, 0.9, False)

render("exp", lp3, sig, cutoff, 0, True)
render("exp", lp3, sig, cutoff, 0.5, True)
render("exp", lp3, sig, cutoff, 0.9, True)

cutoff = np.hstack(
    [peak * np.linspace(1, 0, int(duration * samplerate)) for peak in maxCutoff])

render("lin", lp3, sig, cutoff, 0, False)
render("lin", lp3, sig, cutoff, 0.5, False)
render("lin", lp3, sig, cutoff, 0.9, False)

render("lin", lp3, sig, cutoff, 0, True)
render("lin", lp3, sig, cutoff, 0.5, True)
render("lin", lp3, sig, cutoff, 0.9, True)
