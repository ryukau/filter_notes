import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile

from collections import deque

def nextTime(rng, rate):
    """
    Poisson process.
    https://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/

    If it happens once per N [someunit], then rate is set to 1/N.
    """
    return -np.log(1.0 - rng.uniform(0, 1)) / rate

def pulseNoise(rng, rate, length):
    out = np.zeros(length)
    time = 0
    while time < length:
        time += nextTime(rng, rate)

        t1 = int(time)
        t2 = t1 + 1
        amp = rng.uniform(0, 1)
        if t1 < length:
            out[t1] = amp * (t2 - time)
        if t2 < length:
            out[t2] = amp * (time - t1)
    return out

def makeTriangleFir(delay):
    if delay == 0:
        return np.array([1])
    if delay == 1:
        return np.array([0, 1])
    fir = np.interp(
        np.linspace(0, 1, delay + 1),
        [0, 0.5, 1],
        [0, 1, 0],
    )
    return fir / np.sum(fir)

def peakHold(sig, holdTime, neutral=0):
    out = np.empty(len(sig))
    buffer = deque([0 for _ in range(holdTime)])
    hold = deque([])
    for i in range(len(sig)):
        x0 = sig[i]

        # Remove existing hold values which is less than new hold value.
        if len(hold) > 0:
            idx = len(hold) - 1
            while idx >= 0:
                if hold[idx] < x0:
                    hold.pop()
                else:
                    break
                idx -= 1

        hold.append(x0)  # Add new hold value to hold buffer.

        buffer.append(x0)
        delayOut = buffer.popleft()
        if len(hold) > 0 and delayOut == hold[0]:
            hold.popleft()

        out[i] = hold[0] if len(hold) > 0 else neutral

    return out

def renderDoubleAverageFilterTestSignal():
    samplerate = 48000
    hold = 32
    delay = 32

    seed = np.random.default_rng().integers(0, np.iinfo(np.int64).max)
    with open("synth_seed", "w", encoding="utf-8") as fi:
        fi.write(f"{seed}\n")

    rng = np.random.default_rng(seed)
    triangleFir = makeTriangleFir(delay)

    sig = pulseNoise(rng, 1 / 8, 1024)
    holded = peakHold(sig, hold)
    smoothed = signal.lfilter(triangleFir, 1, holded)

    soundfile.write("snd/input.wav", sig, samplerate, subtype="FLOAT")
    soundfile.write("snd/hold.wav", holded, samplerate, subtype="FLOAT")
    soundfile.write("snd/smoothed.wav", smoothed, samplerate, subtype="FLOAT")

renderDoubleAverageFilterTestSignal()
