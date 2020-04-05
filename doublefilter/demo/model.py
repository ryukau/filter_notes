import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as pyplot

class Model:
    def __init__(self, k1, k2):
        self.k1 = k1
        self.k2 = k2

        self.acc1 = 0
        self.vel1 = 0
        self.pos1 = 0

        self.acc2 = 0
        self.vel2 = 0
        self.pos2 = 0

        self.x1 = 0

    def processHP(self, x0):
        self.acc2 = self.k2 * (self.pos1 - self.pos2)
        self.acc1 = -self.k1 * self.pos1 - self.acc2

        self.vel1 += self.acc1
        self.pos1 += self.vel1

        self.vel2 += self.acc2 + x0 - self.x1
        self.pos2 += self.vel2

        self.x1 = x0
        return self.pos2

    def processLP(self, x0):
        self.acc2 = self.k2 * (self.vel1 - self.vel2)
        self.vel2 += self.acc2 + x0 - self.x1
        self.pos2 += self.vel2 * self.k2  #* np.sqrt(self.k1)

        self.acc1 = -self.k1 * self.pos1 - self.acc2
        self.vel1 += self.acc1
        self.pos1 += self.vel1

        self.x1 = x0
        # return self.pos1 * (self.k1 / self.k2)
        return self.pos2

model = Model(0.1, 0.3)

samplerate = 48000
duration = 0.1
frequency = 60

time = np.linspace(0, frequency * duration, int(duration * samplerate))
phase = 2 * np.pi * time
sig = signal.square(phase)
# sig = signal.sawtooth(phase)
# sig = np.sin(phase)

filtered = np.array([model.processLP(x) for x in sig])

pyplot.figure(figsize=(8, 3))
pyplot.plot(time, sig, color="blue", alpha=0.75, lw=1, label="Input")
pyplot.plot(time, filtered, color="red", alpha=0.75, lw=1, label="Filtered")
pyplot.title(f"Filtered Sin Wave ({frequency}Hz)")
pyplot.xlabel("Time")
pyplot.ylabel("Amplitude")
pyplot.grid()
pyplot.legend(loc=1)
pyplot.show()
