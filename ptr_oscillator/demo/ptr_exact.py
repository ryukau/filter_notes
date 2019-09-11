import numpy
import soundfile
from pathlib import Path

# yapf: disable
def PTR1_exact(phi, T, h):
    n = phi / T
    if T<=phi:
        return 2*T*n-T-1
    if 0<=phi and phi<T:
        return -2*h*n+2*T*n-h**2/T+h/T+2*h-T-1
    return 0.0

def PTR2_exact(phi, T, h):
    n = phi / T
    if 2*T<=phi:
        return 2*T*n-2*T-1
    if 0<=phi and phi<T:
        return -h*n**2-(h**2*n)/T+(h*n)/T+2*T*n-h**3/(3*T**2)+h**2/(2*T**2)-h/(6*T**2)+2*h-2*T-1
    if T<=phi and phi<2*T:
        return h*n**2+(h**2*n)/T-(h*n)/T-4*h*n+2*T*n+h**3/(3*T**2)-(2*h**2)/T-h**2/(2*T**2)+(2*h)/T+h/(6*T**2)+4*h-2*T-1
    return 0.0

def PTR3_exact(phi, T, h):
    n = phi / T
    if 3*T<=phi:
        return 2*T*n-3*T-1
    if 0<=phi and phi<T:
        return -(h*n**3)/3-(h**2*n**2)/(2*T)+(h*n**2)/(2*T)-(h**3*n)/(3*T**2)+(h**2*n)/(2*T**2)-(h*n)/(6*T**2)+2*T*n-h**4/(12*T**3)+h**3/(6*T**3)-h**2/(12*T**3)+2*h-3*T-1
    if T<=phi and phi<2*T:
        return (2*h*n**3)/3+(h**2*n**2)/T-(h*n**2)/T-3*h*n**2+(2*h**3*n)/(3*T**2)-(3*h**2*n)/T-(h**2*n)/T**2+(3*h*n)/T+(h*n)/(3*T**2)+3*h*n+2*T*n+h**4/(6*T**3)-h**3/T**2-h**3/(3*T**3)+(3*h**2)/(2*T)+(3*h**2)/(2*T**2)+h**2/(6*T**3)-(3*h)/(2*T)-h/(2*T**2)+h-3*T-1
    if 2*T<=phi and phi<3*T:
        return -(h*n**3)/3-(h**2*n**2)/(2*T)+(h*n**2)/(2*T)+3*h*n**2-(h**3*n)/(3*T**2)+(3*h**2*n)/T+(h**2*n)/(2*T**2)-(3*h*n)/T-(h*n)/(6*T**2)-9*h*n+2*T*n-h**4/(12*T**3)+h**3/T**2+h**3/(6*T**3)-(9*h**2)/(2*T)-(3*h**2)/(2*T**2)-h**2/(12*T**3)+(9*h)/(2*T)+h/(2*T**2)+9*h-3*T-1
    return 0.0
# yapf: enable

class PTROscillator:
    def __init__(self, ptr_func, samplerate, frequency):
        self.ptr_func = ptr_func
        self.samplerate = samplerate
        self.setFrequency(frequency)
        self.phi = 0.0
        self.h = 1.0

    def setFrequency(self, frequency):
        self.T = frequency / self.samplerate

    def setPhase(self, phi, T):
        self.phi = phi
        ratio = self.T / T
        self.h = ratio - numpy.floor(ratio)

    def process(self):
        self.phi += self.T
        if self.phi >= 1.0:
            self.h = 1.0
            self.phi -= 1.0
        return self.ptr_func(self.phi, self.T, self.h)

class Pulsar:
    def __init__(self, samplerate, frequency):
        self.samplerate = samplerate
        self.T = frequency / samplerate
        self.phi = 0.0

    def setFrequency(self, frequency):
        self.T = frequency / self.samplerate

    def process(self):
        self.phi += self.T
        if self.phi >= 1.0:
            self.phi -= 1.0
            return 1.0
        return 0.0

duration = 4.0
samplerate = 44100
frequency = 1000
sync_ratio = 3 / 2

osc = PTROscillator(PTR3_exact, samplerate, sync_ratio * frequency)
pulsar = Pulsar(samplerate, frequency)

wav = numpy.empty(int(duration * samplerate))
freq = numpy.geomspace(frequency / 8, frequency, len(wav))

for i in range(len(wav)):
    if pulsar.process() > 0.5:
        osc.setPhase(pulsar.phi, pulsar.T)
    pulsar.setFrequency(freq[i])
    wav[i] = osc.process()

snd_dir = Path("snd")
prefix = f"sync_exact{sync_ratio}"
if not snd_dir.exists():
    snd_dir.mkdir(parents=True)

soundfile.write(str(snd_dir / f"{prefix}.wav"), wav, samplerate, subtype="FLOAT")
