import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile

from pathlib import Path
from multiprocessing import Pool
from collections import deque

def plotSignal(signals, labels, title, legendLocation=0):
    plt.figure(figsize=(6, 3.5), tight_layout=True)
    plt.title(title)
    cmap = plt.get_cmap("viridis")
    for idx, (sig, label) in enumerate(zip(signals, labels)):
        plt.plot(sig, lw=1, color=cmap(idx / len(signals)), label=label)
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.legend(loc=legendLocation)
    plt.grid()

    filename = "img/" + title.lower().replace(" ", "_") + ".svg"
    plt.savefig(filename)

def peakHoldForwardSimple(sig, init=0):
    out = np.empty(len(sig))
    hold = init
    counter = 0
    for i in range(len(sig)):
        if sig[i] > hold:
            hold = sig[i]
        out[i] = hold
    return out

def peakHold(sig, holdtime, indices, reset=1):
    out = np.empty(len(sig))
    hold = reset
    counter = 0
    for i in indices:
        if counter > 0:
            counter -= 1
        else:
            hold = reset
        if hold <= sig[i]:  # カウンタをリセットするため必ず <= を使う。
            hold = sig[i]
            counter = holdtime
        out[i] = hold
    return out

def peakHoldForward(sig, holdtime, reset=1):
    return peakHold(sig, holdtime, range(len(sig)), reset)

def peakHoldBackward(sig, holdtime, reset=1):
    return peakHold(sig, holdtime, reversed(range(len(sig))), reset)

def plotForwardHoldExample():
    rng = np.random.default_rng(0)

    sig = rng.normal(0, 1, 256)
    absed = np.abs(sig)
    holdFSimple = peakHoldForwardSimple(absed)

    decaySig = sig * np.geomspace(1, 1e-2, len(sig))
    decayAbsed = np.abs(decaySig)
    decayHoldFSimple = peakHoldForwardSimple(decayAbsed)
    decayHoldF = peakHoldForward(decayAbsed, 32, 0)

    plotSignal([sig, absed], ["Input", "Abs."], "Input Signal")
    plotSignal([absed, holdFSimple], ["Abs.", "Simple Hold"], "Simple Forward Peak Hold",
               4)

    plotSignal([decaySig, decayAbsed], ["Input", "Abs."],
               "Signal with Decaying Amplitude")
    plotSignal([decayAbsed, decayHoldFSimple], ["Abs.", "Simple Hold"],
               "Simple Forward Peak Hold of Decaying Signal")
    plotSignal([decayAbsed, decayHoldF], ["Abs.", "Forward Hold"],
               "Forward Peak Hold of Decaying Signal")
    plt.show()

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

def plotSmoothingExample():
    rng = np.random.default_rng(3)

    delay = 32
    triangleFir = makeTriangleFir(delay)

    sig = pulseNoise(rng, 1 / 8, 256)
    holdF = peakHoldForward(sig, 32, 0)
    smoothed = signal.lfilter(triangleFir, 1, holdF)

    delayedSig = np.hstack((np.zeros(delay), sig))
    delayedHoldF = np.hstack((np.zeros(delay), holdF))

    plotSignal([sig, holdF, smoothed], ["Input", "Forward Hold", "Smoothed"],
               "Smoothed Hold")
    plotSignal([delayedSig, smoothed], ["Delayed Signal", "Smoothed"],
               "Smoothed Hold and Delayed Signal")
    # plotSignal([delayedSig, delayedHoldF, smoothed],
    #            ["Delayed Signal", "Delayed Hold", "Smoothed"],
    #            "Scaling-up of An protruding Peak")
    plt.show()

def idealPeakHoldNaive(sig, holdtime: int):
    out = np.empty(len(sig))
    buffer = deque([0 for _ in range(holdtime)])
    counter = 0
    for i in range(len(sig)):
        buffer.append(abs(sig[i]))
        buffer.popleft()
        out[i] = max(buffer)
    return out

def idealPeakHoldFast(sig, holdTime, neutral=0):
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

def plotIdealHoldExample():
    rng = np.random.default_rng(3)

    hold = 32
    delay = 32
    triangleFir = makeTriangleFir(delay)

    sig = pulseNoise(rng, 1 / 8, 256)
    naive = idealPeakHoldNaive(sig, hold)
    fast = idealPeakHoldFast(sig, hold)

    smoothedNaive = signal.lfilter(triangleFir, 1, naive)
    smoothedFast = signal.lfilter(triangleFir, 1, fast)

    delayedSig = np.hstack((np.zeros(hold), sig[:-hold]))

    plotSignal([sig, naive, fast], ["Input", "Naive", "Fast"], "Ideal Hold")

    plotSignal([delayedSig, smoothedNaive, smoothedFast],
               ["Delayed Input", "Naive", "Fast"], "Smoothed Ideal Hold")

    plt.show()

def verifyHoldEnvelope(sig, env, name):
    condition = (sig - env) > 8 * np.finfo(np.float64).eps
    if np.any(condition):
        idx = np.where(condition)
        diff = np.take(sig, idx) - np.take(env, idx)
        print(f"Test failed: {name}", diff, sep="\n")

def jobTestHoldEnvelope(sig, ident):
    hold = 32
    delay = 32
    triangleFir = makeTriangleFir(delay)

    naive = idealPeakHoldNaive(sig, hold)
    fast = idealPeakHoldFast(sig, hold)

    smoothedNaive = signal.lfilter(triangleFir, 1, naive)
    smoothedFast = signal.lfilter(triangleFir, 1, fast)

    delayed = np.hstack((np.zeros(hold), sig[:-hold]))

    verifyHoldEnvelope(delayed, smoothedNaive, f"Naive, {ident}, {delay}")
    verifyHoldEnvelope(delayed, smoothedFast, f"Fast, {ident}, {delay}")

def jobTestIdealHoldRandom(seed):
    rng = np.random.default_rng(seed)
    sig = pulseNoise(rng, 1 / 8, 48000)
    jobTestHoldEnvelope(sig, seed)

def testIdealHoldRandom(nTest=1024):
    with Pool() as pool:
        seeds = [i for i in range(nTest)]
        idx = 0
        for result in pool.imap_unordered(jobTestIdealHoldRandom, seeds):
            print(f"\r{idx}", end="")
            idx += 1

def jobTestIdealHoldFile(path):
    if not path.is_file():
        return None

    try:
        data, samplerate = soundfile.read(str(path), always_2d=True)
    except:
        print(f"Skipping: {path}")
        return None

    print(path)
    sig = np.abs(data.T[0])
    jobTestHoldEnvelope(sig, str(path))

def testIdealHoldFile(dataDir):
    if not dataDir.exists():
        print(f"Data directory doesn't exists: {dataDir}")
        return

    paths = [p for p in dataDir.glob("**/*")]
    with Pool() as pool:
        for _ in pool.imap_unordered(jobTestIdealHoldFile, paths):
            pass

if __name__ == "__main__":
    plotForwardHoldExample()
    plotSmoothingExample()
    plotIdealHoldExample()
    # testIdealHoldRandom()
    # testIdealHoldFile(Path("../data"))
