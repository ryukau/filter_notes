import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as pyplot

class Model:
    def __init__(self, damping, spring_k):
        self.c = damping
        self.k = spring_k

        self.acc = 0
        self.vel = 0
        self.pos = 0
        self.x1 = 0

    def process(self, x0):
        self.acc = self.c * self.vel + self.k * self.acc
        self.vel -= self.acc + x0 - self.x1
        self.pos -= self.c / (1 - self.k) * self.vel

        self.x1 = x0
        return self.pos

model = Model(1e-1, 0.9)

samplerate = 48000
duration = 0.1
frequency = 60

time = np.linspace(0, frequency * duration, int(duration * samplerate))
phase = 2 * np.pi * time
sig = signal.square(phase)
# sig = np.sin(phase)

filtered = np.array([model.process(x) for x in sig])

pyplot.figure(figsize=(8, 3))
pyplot.plot(time, sig, color="blue", alpha=0.75, lw=1, label="Input")
pyplot.plot(time, filtered, color="red", alpha=0.75, lw=1, label="Filtered")
pyplot.title(f"Filtered Signal ({frequency}Hz)")
pyplot.xlabel("Time")
pyplot.ylabel("Amplitude")
pyplot.grid()
pyplot.legend(loc=1)
pyplot.show()
