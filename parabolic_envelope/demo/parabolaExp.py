import numpy
import math
import matplotlib.pyplot as pyplot

class AccelEnvelope:
    def __init__(
        self,
        samplerate: float,
        lengthA: float,  # In seconds.
        betaA: float,  # In (0, 1).
        decay: float,  # In seconds.
    ):
        # Palabolic attack.
        self.counter: int = -1
        self.y: float = 0
        self.vy: float = 0

        fs2 = samplerate * samplerate

        self.a_A: float = 2 / (betaA * lengthA * lengthA) / fs2
        self.b_A: float = 2 / ((1 - betaA) * lengthA * lengthA) / fs2
        self.t_A: int = int(lengthA * betaA * samplerate)
        self.t_p: int = int(lengthA * samplerate)

        # Exponential decay.
        self.threshold: float = 1e-5
        self.γ_D: float = numpy.power(self.threshold, 1 / (decay * samplerate))
        self.valueD: float = 1

        # For normalization.
        self.gain = 1 / self.getPeak()
        # self.gain = 1 / self.getPeakAlt(samplerate, decay)  # Same result.

    def getPeak(self):
        self.v_A = 0

        log_γ = numpy.log(self.γ_D)
        if self.t_A * log_γ <= -2:
            self.t_k = -2 / numpy.log(self.γ_D)  # Assign to self.t_k for debug.
            return self.γ_D**self.t_k * self.a_A * self.t_k * self.t_k / 2
        print("here")
        temp = self.γ_D**self.t_A * self.a_A
        v_A = temp * self.t_A * (self.t_A / 2 * log_γ + 1)
        h_A = temp * self.t_A * self.t_A / 2
        self.v_A = v_A

        sqrt = numpy.sqrt((v_A * v_A + 2 * self.b_A * h_A) * log_γ * log_γ +
                          self.b_A * self.b_A)
        t_k = (sqrt + v_A * log_γ - self.b_A) / (self.b_A * log_γ)
        self.t_k = t_k + self.t_A

        print(
            sqrt,
            v_A * v_A,
            2 * self.b_A * h_A,
            (v_A * v_A + 2 * self.b_A * h_A) * log_γ * log_γ,
            self.b_A * self.b_A,
        )

        return self.γ_D**t_k * (-self.b_A * t_k * t_k / 2 + v_A * t_k + h_A)

    def getPeakAlt(self, samplerate, decay):
        """Alternative computation."""
        γ = numpy.log(self.threshold) / (decay * samplerate)
        if γ * self.t_A <= -2:
            self.t_k = -2 / γ
            return numpy.exp(γ * self.t_k) * self.a_A * self.t_k * self.t_k / 2
        temp = numpy.exp(γ * self.t_A) * self.a_A
        v_A = temp * self.t_A * (γ * self.t_A / 2 + 1)
        h_A = temp * self.t_A * self.t_A / 2

        sqrt = numpy.sqrt((v_A * v_A + 2 * self.b_A * h_A) * γ * γ + self.b_A * self.b_A)
        self.t_k = (sqrt + v_A * γ - self.b_A) / (self.b_A * γ)

        return numpy.exp(
            γ * self.t_k) * (-self.b_A * self.t_k * self.t_k / 2 + v_A * self.t_k + h_A)

    def process(self):
        self.counter += 1
        if self.counter < self.t_A:
            self.vy += self.a_A
            self.y += self.vy
        elif self.counter < self.t_p:
            self.vy -= self.b_A
            self.y += self.vy
        self.valueD *= self.γ_D
        return self.gain * self.y * self.valueD

def plotEnvelope():
    samplerate = 48000
    lengthA = 2
    betaA = 0.2
    decay = 4

    envelope = AccelEnvelope(samplerate, lengthA, betaA, decay)

    duration = decay
    time = numpy.linspace(0, duration, duration * samplerate)
    data = [envelope.process() for _ in range(len(time))]

    pyplot.figure(figsize=(8, 3))
    pyplot.plot(time, data, color="black", label="Envelope")
    if envelope.t_A / samplerate < duration:
        pyplot.axvline(
            envelope.t_A / samplerate, color="#ff0000", lw=1, ls="--", label="t_A")
        pyplot.scatter(envelope.t_A / samplerate, envelope.v_A * envelope.gain)
    pyplot.axvline(envelope.t_k / samplerate, color="#0044ff", lw=1, ls="--", label="t_k")
    if envelope.t_p / samplerate < duration:
        pyplot.axvline(
            envelope.t_p / samplerate, color="#ff8800", lw=1, ls="--", label="t_p")
    pyplot.title(
        f"Parabolic Attack * Exponential Decay (l_A={lengthA}, β={betaA:.1f}, t_d={decay})"
    )
    pyplot.xlabel("Time [s]")
    pyplot.ylabel("Amplitude")
    pyplot.grid()
    pyplot.legend()
    pyplot.tight_layout()
    pyplot.show()

def plotPeak():
    samplerate = 48000
    betaA = numpy.linspace(0.01, 0.99, 11)
    lengthA = numpy.linspace(0.01, 0.99, 11)
    decay = 1

    cmap = pyplot.get_cmap("plasma")
    pyplot.figure(figsize=(8, 5))
    for idx, beta in enumerate(betaA):
        peaks = []
        for length in lengthA:
            envelope = AccelEnvelope(samplerate, length, beta, decay)
            data = [envelope.process() for _ in range(int(length * samplerate))]
            peaks.append(numpy.max(data))

        pyplot.plot(
            lengthA,
            peaks,
            alpha=0.75,
            color=cmap(idx / len(betaA)),
            label=f"β={beta:.3f}")
    pyplot.title("Real Peak")
    pyplot.xlabel("AttackTime [s]")
    pyplot.ylabel("Peak")
    pyplot.grid()
    pyplot.legend(loc="lower left", ncol=2)
    pyplot.tight_layout()
    pyplot.show()

plotEnvelope()
# plotPeak()
