import numpy
import matplotlib.pyplot as pyplot

class AccelEnvelope:
    def __init__(
        self,
        samplerate: float,
        lengthA: float,  # In seconds.
        betaA: float,  # In (0, 1).
        lengthR: float,  # In seconds.
        betaR: float,  # In (0, 1).
    ):
        self.counter: int = -1
        self.y: float = 0
        self.vy: float = 0

        fs2 = samplerate * samplerate

        self.a_A: float = 2 / (betaA * lengthA * lengthA) / fs2
        self.b_A: float = 2 / ((1 - betaA) * lengthA * lengthA) / fs2
        self.n_A: int = int(lengthA * betaA * samplerate)
        self.n_p: int = int(lengthA * samplerate)

        self.a_R: float = 2 / (betaR * lengthR * lengthR) / fs2
        self.b_R: float = 2 / ((1 - betaR) * lengthR * lengthR) / fs2
        self.n_R: int = self.n_p + int(lengthR * betaR * samplerate)
        self.n_E: int = self.n_p + int(lengthR * samplerate)

    def process(self):
        self.counter += 1
        if self.counter < self.n_A:
            self.vy += self.a_A
        elif self.counter < self.n_p:
            self.vy -= self.b_A
        elif self.counter == self.n_p:  # Just in case.
            self.vy = 0
            return self.y
        elif self.counter < self.n_R:
            self.vy -= self.a_R
        elif self.counter < self.n_E and self.y > 1e-5:
            self.vy += self.b_R
        else:
            return 0
        self.y += self.vy
        return self.y

samplerate = 48000
lengthA = 2
betaA = 0.2
lengthR = 3
betaR = 0.2

envelope = AccelEnvelope(samplerate, lengthA, betaA, lengthR, betaR)

duration = 6
time = numpy.linspace(0, duration, duration * samplerate)
data = [envelope.process() for _ in range(len(time))]

pyplot.figure(figsize=(8, 3))
pyplot.plot(time, data, color="black", label="Envelope")
pyplot.axvline(envelope.n_A / samplerate, color="#ff0000", lw=1, ls="--", label="t_A")
pyplot.axvline(envelope.n_p / samplerate, color="#ff4400", lw=1, ls="--", label="t_p")
pyplot.axvline(envelope.n_R / samplerate, color="#ff8800", lw=1, ls="--", label="t_R")
pyplot.axvline(envelope.n_E / samplerate, color="#ffbb00", lw=1, ls="--", label="t_E")
pyplot.title("Parabolic Envelope")
pyplot.xlabel("Time [s]")
pyplot.ylabel("Amplitude")
pyplot.grid()
pyplot.legend()
pyplot.tight_layout()
pyplot.show()
