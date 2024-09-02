import matplotlib.pyplot as pyplot
import numpy
import soundfile
import numpy as np
from scipy.special import lambertw
from dataclasses import dataclass


@dataclass
class EnvelopeResult:
    raw: np.array
    normalized: np.array
    peakTime: float


def envelope1(
    attackSeconds: float,
    decaySeconds: float,
    time: np.array,
    eps=np.finfo(np.float64).eps,
):
    """`attackSeconds > 0` and `decaySeconds > 0`."""

    def env(a, d, t):
        return (1 - np.exp(a * t)) * np.exp(d * t)

    a = np.log(eps) / attackSeconds
    d = np.log(eps) / decaySeconds
    peakTime = -np.log1p(a / d) / a
    gain = 1 / env(a, d, peakTime)

    # print(f"E1: t_p={peakTime}, a={a}, d={d}")  # debug

    raw = env(a, d, time)
    return EnvelopeResult(raw, raw * gain, peakTime)


def envelope2(
    attackSeconds: float,
    decaySeconds: float,
    time: np.array,
    eps=np.finfo(np.float64).eps,
):
    """`attackSeconds > 0` and `decaySeconds > 0`."""

    def env(a, d, t):
        return (1 - a**t) * d**t

    a = np.power(eps, 1 / attackSeconds)
    d = np.power(eps, 1 / decaySeconds)
    peakTime = -np.log1p(a / d) / a

    log_a = np.log(a)
    log_d = np.log(d)
    peakTime = np.log(log_d / (log_d + log_a)) / log_a
    gain = 1 / env(a, d, peakTime)

    # print(f"E2: t_p={peakTime}, a={a}, d={d}")  # debug

    raw = env(a, d, time)
    return EnvelopeResult(raw, raw * gain, peakTime)


def envelope3(
    peakSeconds: float,
    releaseSeconds: float,
    time: np.array,
    eps=np.finfo(np.float64).eps,
):
    """
    `peakSeconds > 0`.
    `releaseSeconds` adds to the release time calculated from `peakSeconds`.
    """

    def env(a, d, t):
        return (1 - np.exp(a * t)) * np.exp(d * t)

    D = releaseSeconds - np.log(eps) * peakSeconds
    d = np.log(eps) / D
    x = d * peakSeconds
    a = (lambertw(x * np.exp(x), -1) / peakSeconds - d).real
    gain = 1 / env(a, d, peakSeconds)

    # # debug
    # print(f"E3: t_p={peakSeconds}, a={a}, d={d}")
    # print(f"A={-np.log(eps)/np.log(-a)}, D={D}")
    # print(np.exp(a))
    # print(np.exp(d))

    raw = env(a, d, time)
    return EnvelopeResult(raw, raw * gain, peakSeconds)


def compareExpAD(name, targetFunc, arg1, arg2, epsilon):
    data, samplerate = soundfile.read(f"snd/{name}.wav", always_2d=True)
    data = data.T[0]

    time = np.linspace(0, len(data) / samplerate, len(data), endpoint=False)
    target = targetFunc(arg1, arg2, time, epsilon)

    pyplot.figure(figsize=(6, 3))
    pyplot.title(f"Test of {name}")
    pyplot.axvline(target.peakTime, color="black", alpha=0.5, lw=1, ls=":", label="t_p")
    pyplot.plot(time, target.normalized, color="blue", alpha=0.75, label="Target")
    pyplot.plot(time, data, color="red", alpha=0.75, label="C++")
    pyplot.xlabel("Time [s]")
    pyplot.ylabel("Amplitude")
    pyplot.legend()
    pyplot.grid()
    pyplot.tight_layout()
    pyplot.show()


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
    phase = numpy.linspace(
        0, 2 * numpy.pi * frequency * len(data) / samplerate, len(data)
    )
    tone = numpy.sin(phase)
    soundfile.write(f"snd/tone_{name}.wav", data * tone, samplerate, subtype="FLOAT")


process("ADSR")
process("ReleaseWhileAttack")
process("ReleaseWhileDecay")
process("TriggerWhileRelease")
process("ChangeSustain")

process("P_ADSR")
process("P_ReleaseWhileAttack")
process("P_ReleaseWhileDecay")
process("P_TriggerWhileRelease")
process("P_ChangeSustain")

process("ExpAD1")
process("ExpAD2")
process("ExpAD3")

compareExpAD("ExpAD1", envelope1, 1, 2, 1e-5)
compareExpAD("ExpAD2", envelope2, 1, 2, 1e-5)
compareExpAD("ExpAD3", envelope3, 0.1, 2, np.finfo(np.float64).eps)
