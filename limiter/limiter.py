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

def doubleMovingAverageFilter(sig, delay):
    hd = delay // 2

    delay1 = deque([0 for _ in range(hd)])
    delay2 = deque([0 for _ in range(hd)])
    sum1 = 0
    sum2 = 0
    buf = 0

    out = np.empty_like(sig)
    for i in range(len(sig)):
        delay1.append(buf)
        sum1 += buf - delay1.popleft()
        out1 = sum1 / hd

        delay2.append(out1)
        sum2 += out1 - delay2.popleft()
        out[i] = sum2 / hd

        buf = sig[i]
    return out

def plotDMAFilter():
    delay = 32

    sig = np.zeros(delay * 2 + 1)
    sig[0] = 1

    triangleFir = makeTriangleFir(delay)
    smoothed = signal.lfilter(triangleFir, 1, sig)

    dma = doubleMovingAverageFilter(sig, delay)

    plt.plot(smoothed, label="Ref.")
    plt.plot(dma, label="Fast")
    plt.grid()
    plt.legend()
    plt.show()

def limiter():
    samplerate = 48000
    rng = np.random.default_rng(3)

    holdtime = 32
    delay = 32
    th = 0.05
    threshold = th - 1e-15
    triangleFir = makeTriangleFir(delay)

    sig = pulseNoise(rng, 1 / 8, 256)
    holded = peakHold(np.abs(sig), holdtime)
    candidate = np.where(holded > threshold, threshold / holded, 1)
    gain = candidate  # TODO Add release.
    delayed = np.hstack((np.zeros(delay), sig[:-delay]))

    smoothedA = signal.lfilter(triangleFir, 1, gain)
    outA = smoothedA * delayed

    smoothedB = doubleMovingAverageFilter(gain, delay)
    outB = smoothedB * delayed

    if np.any(outA > th):
        print("Test A failed.")
    if np.any(outB > th):
        print("Test B failed.")

    # plt.plot(outA)
    # plt.plot(outB)
    plt.plot(smoothedA, label="A")
    plt.plot(smoothedB, label="B")
    plt.plot(gain, label="gain")
    plt.grid()
    plt.legend()
    plt.show()

def plotHardClipCharacteristic():
    tick = np.linspace(0, 1.2, 128)
    character = np.where(tick > 1, 1, tick)

    plt.figure(figsize=(3.5, 3.5), tight_layout=True)
    plt.title("Hard Clip Curve")
    plt.plot(tick, character)
    plt.xlabel("Input Amplitude")
    plt.ylabel("Output Amplitude")
    plt.grid()
    plt.show()

def softclip(x0, threshold=1, ratio=0.9):
    h = threshold
    absed = np.abs(x0)

    rh = h * ratio
    if absed <= rh:
        return x0

    xc = 2 * h - rh
    if absed >= xc:
        return h

    return np.sign(x0) * (h + (xc - absed) * (xc - absed) * 0.25 / (rh - h))

def testSoftClip():
    sig = np.linspace(0, 1.5, 1024)
    out = np.empty_like(sig)
    for i in range(len(sig)):
        out[i] = softclip(sig[i])

    plt.plot(out)
    plt.grid()
    plt.show()

testSoftClip()
