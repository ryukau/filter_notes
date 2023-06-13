import math
import numpy as np
import matplotlib.pyplot as plt
import time
import json
from multiprocessing import Pool
from scipy.stats import linregress


class FeedbackEMALowpass:
    def __init__(self, sampleRate, cutoffHz, resonance):
        self.reset(resonance)

        cutoff = np.clip(cutoffHz / sampleRate, 0, 0.5)

        y = 1 - math.cos(math.tau * cutoff)
        self.c1 = math.sqrt((y + 2) * y) - y

        t = math.tan(math.pi * cutoff)
        self.c2 = (t - 1) / (t + 1)

    def reset(self, resonance):
        self.u1 = 0
        self.v1 = 0
        self.u2 = 0

        self.q = resonance

    def process(self, x0):
        self.v1 = self.c2 * (self.u1 - self.v1) + self.u2
        self.u2 = self.u1
        self.u1 += self.c1 * (x0 - self.u1) - self.q * self.v1
        return self.u1


def plotIR(sampleRate, cutoffHz, resonance=0.5):
    # nSample = int(8 * sampleRate / cutoffHz)
    nSample = sampleRate

    sig = np.zeros(nSample)

    flt = FeedbackEMALowpass(sampleRate, cutoffHz, resonance)
    sig[0] = flt.process(1)
    for i in range(1, nSample):
        sig[i] = flt.process(0)

    plt.plot(np.abs(sig))
    plt.axhline(1, color="red", lw=1)
    plt.yscale("log")
    plt.grid()
    plt.show()


def find1(sampleRate: int, cutoffHz):
    nSample = sampleRate

    # Binary search.
    high = 2
    low = 0
    mid = 1
    iteration = 0
    flt = FeedbackEMALowpass(sampleRate, cutoffHz, mid)
    while True:
        flt.reset(mid)
        flt.process(1)

        index = 0
        while index < nSample:
            sig = flt.process(0)
            if sig > 1 or sig < -1:
                high = mid
                break
            index += 1
        if index >= nSample:
            low = mid
        mid = (high + low) / 2
        print(low, mid, high)

        if np.abs(high - low) <= np.finfo(np.float64).eps:
            break
        iteration += 1

    print(f"iteration: {iteration}, cutoff: {cutoffHz/sampleRate}, max Q: {low}")
    plotIR(sampleRate, cutoffHz, low)


def isConversing(sampleRate, cutoffHz, resonance):
    flt = FeedbackEMALowpass(sampleRate, cutoffHz, resonance)
    y1 = 0
    y0 = flt.process(1)
    while y1 < y0:
        y1 = y0
        y0 = flt.process(0)
    peak1 = y1
    while y1 > y0:
        y1 = y0
        y0 = flt.process(0)
    peak2 = abs(y1)

    if math.isfinite(peak1) and math.isfinite(peak2):
        return peak1 >= peak2
    return False


def find2(sampleRate, cutoffHz):
    # Binary search.
    high = 1
    low = 0
    mid = 0.5
    iteration = 0
    while True:
        if isConversing(sampleRate, cutoffHz, mid):
            low = mid
        else:
            high = mid
        mid = (low + high) / 2

        if np.abs(high - low) < np.finfo(np.float64).eps:
            break
        iteration += 1

    print(f"iteration: {iteration}, cutoff: {cutoffHz/sampleRate}, max Q: {mid}")
    # plotIR(sampleRate, cutoffHz, mid)
    return low


def find3(sampleRate: int, cutoffHz):
    nSample = sampleRate if cutoffHz / sampleRate >= 0.02 else 4 * sampleRate

    # Binary search.
    high = 2
    low = 0
    mid = 1
    iteration = 0
    flt = FeedbackEMALowpass(sampleRate, cutoffHz, mid)
    time = np.arange(nSample)
    slope = None
    while True:
        flt.reset(mid)

        sig = np.zeros(nSample)
        flt.process(1)
        isDiversing = False
        for i in range(nSample):
            sig[i] = flt.process(0)
            if sig[i] > 1 or sig[i] < -1:
                isDiversing = True
                break

        if isDiversing:
            high = mid
        else:
            slope = linregress(time, np.abs(sig)).slope
            if slope < 0:
                low = mid
            elif slope > 0:
                high = mid
            else:
                print("here")
                break
        mid = (low + high) / 2

        if np.abs(high - low) <= np.finfo(np.float64).eps:
            break
        iteration += 1

    print(f"iteration: {iteration}, cutoff: {cutoffHz}, max Q: {low}")
    # plotIR(sampleRate, cutoffHz, low)
    return low


def job(sampleRate, cutoffHz):
    return (cutoffHz / sampleRate, find3(sampleRate, cutoffHz))


def getCurveLinear(sampleRate):
    data = {"normalizedFreq": [], "resonance": []}

    args = [(sampleRate, cutoffHz) for cutoffHz in range(1, sampleRate // 2)]
    with Pool() as pool:
        result = pool.starmap(job, args)

    result = list(map(list, zip(*result)))
    data = {"normalizedFreq": [0] + result[0], "resonance": [0] + result[1]}
    with open("table_linear.json", "w", encoding="utf-8") as fi:
        json.dump(data, fi)


def getCurveLog(sampleRate):
    data = {"normalizedFreq": [], "resonance": []}

    midinote = np.linspace(0, 136, sampleRate // 2)
    hertz = 440 * np.exp2((midinote - 69) / 12)
    args = [(sampleRate, cutoffHz) for cutoffHz in hertz]
    with Pool() as pool:
        result = pool.starmap(job, args)

    result = list(map(list, zip(*result)))
    data = {"normalizedFreq": result[0], "resonance": result[1]}
    with open("table_log.json", "w", encoding="utf-8") as fi:
        json.dump(data, fi)


if __name__ == "__main__":
    # plotIR(48000, 20000, 1.0604970675703127)

    # start = time.time()
    # find1(48000, 20000)
    # end = time.time()
    # print(f"Elapsed Time [s]: {end - start}")

    # start = time.time()
    # find2(48000, 1000)
    # end = time.time()
    # print(f"Elapsed Time [s]: {end - start}")

    # start = time.time()
    # find3(48000, 15555)
    # end = time.time()
    # print(f"Elapsed Time [s]: {end - start}")

    getCurveLinear(65536)
    getCurveLog(65536)
