import numpy as np
import matplotlib.pyplot as plt
import soundfile
import functools

from pathlib import Path

def sinwave(samplerate=48000, duration=1, frequency=100, gainDb=-6):
    nFrames = int(samplerate * duration)
    phase = np.linspace(0, 2 * np.pi * frequency * duration, nFrames)
    return np.power(10, gainDb / 20) * np.sin(phase)

def envelope(samplerate, duration, cycle=0.1):
    env = np.geomspace(1, 1e-5, int(samplerate * cycle)) - 1e-5
    length = int(samplerate * duration)
    tiled = np.tile(env, length // len(env))
    return np.hstack((
        tiled,
        np.zeros(length - len(tiled)),
    ))

def fade(sig, ratio):
    fadeFrame = int(len(sig) * ratio)
    flatFrame = len(sig) - 2 * fadeFrame
    envelope = np.hstack((
        np.linspace(0, 1, fadeFrame),
        np.ones(flatFrame),
        np.linspace(1, 0, fadeFrame),
    ))
    return sig * envelope

def panMonoLinear(pan):
    """
    pan == 0.0: left
    pan == 0.5: center
    pan == 1.0: right

    return (gainL, gainR)
    """
    return (1 - pan, pan)

def panMonoConstantPower(pan):
    theta = (0.5 * np.pi) * pan
    return (np.cos(theta), np.sin(theta))

def panMonoIntermediatePower(pan):
    """-4.5 dB power panning law"""
    theta = (0.5 * np.pi) * pan
    return (
        np.sqrt((1 - pan) * np.cos(theta)),
        np.sqrt(pan * np.sin(theta)),
    )

def getStereoGain(gainFunc, param):
    LL = gainFunc(pan, param)
    RR = gainFunc(1 - pan, param)
    RL = 1 - LL
    LR = 1 - RR
    return (LL, RL, LR, RR)

def panStereoLinear(pan, param=None):
    def gain(pan, unused):
        return np.where(pan >= 0.5, 1, 0.5 + pan)

    return getStereoGain(gain, param)

def panStereoPartial2ndOrder(pan, param=0.8):
    """`param` is in [0.0, 1.0]."""
    def gain(p, a):
        h = 1 / (a + 0.5)
        return np.where(
            p <= a,
            h * p + 0.5,
            np.where(
                p >= 0.5,
                1,
                h * (p - 0.5)**2 / (2 * a - 1) + 1,
            ),
        )

    return getStereoGain(gain, 0.5 * param)

def panStereoPartialSin(pan, param=0.8):
    """`param` is in [0.0, 1.0]."""
    def gain(p, a):
        h = 2 / (2 * a + 1)
        return np.where(
            p <= a,
            h * p + 0.5,
            np.where(
                p >= 0.5,
                1,
                0.5 * h * ((0.5 - a) / np.pi * np.sin(
                    (p - a) / (0.5 - a) * np.pi) + p - 0.5) + 1,
            ),
        )

    return getStereoGain(gain, 0.5 * param)

def panStereoCircle(pan, param=1):
    """`param` is in (0.0, 1.0]."""
    def gain(p, x):
        y = np.sqrt(1 - x * x)
        return np.where(
            p <= 0.5,
            0.5 + 0.5 * np.sqrt(1 - (x - 2 * x * p)**2) - y / (1 - y),
            1,
        )

    return getStereoGain(gain, param)

def panStereoPoly(pan, param=0.25):
    """`param` is in [0.0, 1.0]."""
    def gain(p, n):
        return np.where(
            p <= 0.5,
            1 - 0.5 * np.power(1 - 2 * p, n),
            1,
        )

    order = (1 + 4 * param)
    return getStereoGain(gain, order)

def panStereoSin(pan, param=None):
    def gain(p, unused):
        return np.where(
            p <= 0.5,
            0.5 + 0.5 * np.sin(np.pi * p),
            1,
        )

    return getStereoGain(gain, param)

def panStereoSCurve(pan, param=0):
    def gain(p, n):
        return np.where(
            p <= 0.5,
            0.75 - 0.25 * np.cos(2 * np.pi * p)**n,
            1,
        )

    order = 1 + 3 * param
    return getStereoGain(gain, order)

def panStereoSoftplus(pan, param=0.5):
    """`param` is in [0.0, 1.0]."""
    def softplus(x):
        return np.log(1 + np.exp(x))

    def gain(p, k):
        v_min = softplus(-0.5 * k)
        v_max = softplus(0.5 * k)
        return 0.5 + 0.5 * (v_max - softplus(k * (0.5 - p))) / (v_max - v_min)

    saturation = np.power(10, 2 * param)
    return getStereoGain(gain, saturation)

def panStereoSinc(pan, param=0.5):
    """`param` is in [0.0, 1.0]."""
    def gain(p, k):
        return np.where(
            p <= 0.5,
            0.5 + 0.5 * np.sinc(k * (1 - 2 * p)),
            1,
        )

    zeroCross = 1 + np.floor(16 * param)
    return getStereoGain(gain, zeroCross)

def testMono(panFunc, pan, source, samplerate):
    """source must be 1D array."""
    gainL, gainR = panFunc(pan)

    sigL = gainL * source
    sigR = gainR * source

    ## Gain plot.
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches((10, 2.5))
    fig.set_tight_layout(True)
    ax[0].set_title("Gain")
    ax[0].plot(pan, gainL, label="L")
    ax[0].plot(pan, gainR, label="R")
    ax[0].set_xlabel("Pan")
    ax[0].set_ylabel("Gain")
    ax[0].set_xticks(np.linspace(0, 1, 5))
    ax[0].set_ylim((-0.05, 1.05))
    ax[0].grid()
    ax[0].legend()

    ## L + R plot.
    mergedGainDecibel = 20 * np.log10(gainL + gainR)
    ax[1].set_title("L + R")
    ax[1].plot(pan, mergedGainDecibel, label="L + R")
    ax[1].set_xlabel("Pan")
    ax[1].set_ylabel("Gain [dB]")
    ax[1].set_xticks(np.linspace(0, 1, 5))
    ax[1].set_ylim((-0.25, 3.25))
    ax[1].grid()
    ax[1].legend()

    ## Power plot.
    powerDecibel = 10 * np.log10(gainL**2 + gainR**2)
    ax[2].set_title("Power")
    ax[2].plot(pan, powerDecibel, label="L^2 + R^2")
    ax[2].set_xlabel("Pan")
    ax[2].set_ylabel("Power [dB]")
    ax[2].set_xticks(np.linspace(0, 1, 5))
    # ax[2].set_yticks(np.linspace(0, 1, 5))
    ax[2].set_ylim((-3.25, 0.25))
    ax[2].grid()
    ax[2].legend()
    plt.savefig(f"img/{panFunc.__name__}.svg")
    plt.close("all")

    data = np.array((sigL, sigR)).T
    soundfile.write(f"snd/{panFunc.__name__}.wav", data, samplerate, subtype="FLOAT")

def testStereo(panFuncMono, panFuncStereo, pan, source, samplerate):
    """source must be 2D array."""
    LL, RL, LR, RR = panFuncStereo(pan)
    gainL, gainR = panFuncMono(pan)

    sigL = gainL * (LL * source[0] + RL * source[1])
    sigR = gainR * (LR * source[0] + RR * source[1])

    plt.figure(figsize=(4, 2.5), tight_layout=True)
    plt.plot(pan, LL, label="G_LL")
    plt.plot(pan, RR, label="G_RR")
    plt.xlabel("Pan")
    plt.ylabel("Gain")
    plt.xticks(np.linspace(0, 1, 5))
    plt.yticks(np.linspace(0.5, 1, 5))
    plt.grid()
    plt.legend()
    plt.savefig(f"img/{panFuncStereo.__name__}.svg")
    plt.close("all")

    data = np.array((sigL, sigR)).T
    soundfile.write(
        f"snd/{panFuncStereo.__name__}_{panFuncMono.__name__}.wav",
        data,
        samplerate,
        subtype="FLOAT",
    )

def testSignal(samplerate, duration=4, fadeSeconds=0.05):
    ping = envelope(samplerate, duration) * sinwave(samplerate, duration, 440)
    ping = fade(ping, fadeSeconds / duration)

    tone = sinwave(samplerate, duration, 110)
    tone = fade(tone, fadeSeconds / duration)

    return (ping, tone)

if __name__ == "__main__":
    Path("snd").mkdir(parents=True, exist_ok=True)

    samplerate = 48000

    stereo = testSignal(samplerate, duration=8)
    mono = 0.5 * (stereo[0] + stereo[1])

    pan = np.linspace(0, 1, len(mono))

    testM = functools.partial(testMono, pan=pan, source=mono, samplerate=samplerate)
    testM(panMonoLinear)
    testM(panMonoConstantPower)
    testM(panMonoIntermediatePower)

    for panFuncMono in [panMonoLinear, panMonoConstantPower, panMonoIntermediatePower]:
        testS = functools.partial(
            testStereo,
            panFuncMono=panFuncMono,
            pan=pan,
            source=stereo,
            samplerate=samplerate,
        )
        testS(panFuncStereo=panStereoLinear)
        testS(panFuncStereo=panStereoPartial2ndOrder)
        testS(panFuncStereo=panStereoPartialSin)
        testS(panFuncStereo=panStereoHalfCircle)
        testS(panFuncStereo=panStereoPoly)
        testS(panFuncStereo=panStereoSin)
        testS(panFuncStereo=panStereoSCurve)
        testS(panFuncStereo=panStereoSoftplus)
        testS(panFuncStereo=panStereoSinc)
