import numpy
import matplotlib.pyplot as pyplot

class AccelEnvelope:
    """
    timeR is offset from time at peak (t_p). For example, if t_p = 1 and timeR = 0.5,
    accelaration is change at 1.5 seconds from the start of envelope.
    """
    def __init__(
        self,
        samplerate: float,
        accelA: float,
        brakeA: float,
        timeA: float,
        accelR: float,
        brakeR: float,
        timeR: float,
    ):
        self.counter: int = -1
        self.y: float = 0
        self.vy: float = 0

        self.v_A: float = accelA * timeA
        self.h_A: float = accelA * timeA * timeA / 2
        self.t_p: float = self.v_A / brakeA
        self.h_p: float = -brakeA / 2 * self.t_p * self.t_p + self.v_A * self.t_p + self.h_A

        self.v_R: float = accelR * timeR
        self.h_R: float = accelR * timeR * timeR / 2
        self.t_E: float = self.v_R / brakeR
        self.h_E: float = -brakeR / 2 * self.t_E * self.t_E + self.v_R * self.t_E + self.h_R

        accelA, brakeA = AccelEnvelope.normalize(accelA, brakeA, timeA)
        accelR, brakeR = AccelEnvelope.normalize(accelR, brakeR, timeR)

        # Change time unit from seconds to samples.
        fs2 = samplerate * samplerate

        self.accelA: float = accelA / fs2
        self.brakeA: float = brakeA / fs2
        self.timeA: int = int(timeA * samplerate)
        self.n_p: int = int((self.t_p + timeA) * samplerate)

        self.accelR: float = accelR / fs2
        self.brakeR: float = brakeR / fs2
        self.timeR: int = self.n_p + int(timeR * samplerate)
        self.n_E: int = self.n_p + int((self.t_E + timeR) * samplerate)

    def normalize(accel, brake, time):
        nu = 2 * brake / ((accel * brake + accel * accel) * time * time)
        return (accel * nu, brake * nu)

    def process(self):
        self.counter += 1
        if self.counter < self.timeA:
            self.vy += self.accelA
        elif self.counter < self.n_p:
            self.vy -= self.brakeA
        elif self.counter == self.n_p:  # Just in case.
            self.vy = 0
            return self.y
        elif self.counter < self.timeR:
            self.vy -= self.accelR
        elif self.counter < self.n_E and self.y > 1e-5:
            self.vy += self.brakeR
        else:
            return 0
        self.y += self.vy
        return self.y

samplerate = 48000
accelA = 0.001
brakeA = 0.002
timeA = 0.5
accelR = 0.0003
brakeR = 0.0004
timeR = 1.0

envelope = AccelEnvelope(samplerate, accelA, brakeA, timeA, accelR, brakeR, timeR)

duration = 3
time = numpy.linspace(0, duration, duration * samplerate)
timeP = envelope.n_p / samplerate
timeE = envelope.n_E / samplerate
data = [envelope.process() for _ in range(len(time))]

pyplot.figure(figsize=(8, 3))
pyplot.plot(time, data, color="black", label="Envelope")
pyplot.axvline(timeA, color="#ff0000", lw=1, ls="--", label="t_A")
pyplot.axvline(timeP, color="#ff4400", lw=1, ls="--", label="t_p")
pyplot.axvline(timeP + timeR, color="#ff8800", lw=1, ls="--", label="t_R")
pyplot.axvline(timeE, color="#ffbb00", lw=1, ls="--", label="t_E")
pyplot.title("Parabolic Envelope")
pyplot.xlabel("Time [s]")
pyplot.ylabel("Amplitude")
pyplot.grid()
pyplot.legend()
pyplot.tight_layout()
pyplot.show()
