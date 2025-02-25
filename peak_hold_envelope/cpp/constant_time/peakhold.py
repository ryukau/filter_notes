import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import soundfile
from collections import deque


class PeakHold:
    lowest = np.finfo(np.float64).min

    def __init__(S, maxLength: int):
        S.backIndex: int = 0
        S.middleStart: int = 0
        S.workingIndex: int = 0
        S.middleEnd: int = 0
        S.frontIndex: int = 0

        S.frontMax: float = PeakHold.lowest
        S.workingMax: float = PeakHold.lowest
        S.middleMax: float = PeakHold.lowest

        S.resize(maxLength)

    def size(S):
        return S.frontIndex - S.backIndex

    def resize(S, maxLength: int):
        bufferLength = int(2 ** np.ceil(1 + np.log2(maxLength)))
        S.buffer = np.zeros(bufferLength)
        S.bufferMask = bufferLength - 1

        S.frontIndex = S.backIndex + maxLength
        S.reset()

    def reset(S, fill: float = lowest):
        prevSize: int = S.size()
        S.buffer = np.full_like(S.buffer, fill)

        S.frontMax = S.workingMax = S.middleMax = PeakHold.lowest
        S.middleEnd = S.workingIndex = S.frontIndex = 0
        S.middleStart = S.middleEnd - (prevSize // 2)
        S.backIndex = S.frontIndex - prevSize

    def set(S, newSize: int, preserveCurrentPeak: bool = False):
        while S.size() < newSize:
            backPrev: float = S.buffer[S.backIndex & S.bufferMask]
            S.backIndex -= 1

            index: int = S.backIndex & S.bufferMask
            back: float = S.buffer[index]
            S.buffer[index] = backPrev if preserveCurrentPeak else max([back, backPrev])
        while S.size() > newSize:
            S.pop()

    def push(S, v: float):
        S.buffer[S.frontIndex & S.bufferMask] = v
        S.frontIndex += 1
        S.frontMax = max([S.frontMax, v])

    def pop(S):
        if S.backIndex == S.middleStart:
            # Move along the maximums
            S.workingMax = PeakHold.lowest
            S.middleMax = S.frontMax
            S.frontMax = PeakHold.lowest

            prevFrontLength: int = S.frontIndex - S.middleEnd
            prevMiddleLength: int = S.middleEnd - S.middleStart
            if prevFrontLength <= prevMiddleLength + 1:
                # Swap over simply
                S.middleStart = S.middleEnd
                S.middleEnd = S.frontIndex
                S.workingIndex = S.middleEnd
            else:
                # The front is longer than the middle - only happens if unbalanced
                # We don't move *all* of the front over, keeping half the surplus in the front
                middleLength: int = (S.frontIndex - S.middleStart) // 2
                S.middleStart = S.middleEnd
                S.middleEnd += middleLength

                # Working index is close enough that it will be finished by the time the back is empty
                backLength: int = S.middleStart - S.backIndex
                workingLength: int = min([backLength, S.middleEnd - S.middleStart])
                S.workingIndex = S.middleStart + workingLength

                # Since the front was not completely consumed, we re-calculate the front's maximum
                i: int = S.middleEnd
                while i != S.frontIndex:
                    S.frontMax = max([S.frontMax, S.buffer[i & S.bufferMask]])
                    i += 1

                # The index might not start at the end of the working block - compute the last bit immediately
                i: int = S.middleEnd - 1
                while i != S.workingIndex - 1:
                    S.buffer[i & S.bufferMask] = S.workingMax = max(
                        [S.workingMax, S.buffer[i & S.bufferMask]]
                    )
                    i -= 1

            # Is the new back (previous middle) empty? Only happens if unbalanced
            if S.backIndex == S.middleStart:
                # swap over again (front's empty, no change)
                S.workingMax = PeakHold.lowest
                S.middleMax = S.frontMax
                S.frontMax = PeakHold.lowest
                S.middleStart = S.workingIndex = S.middleEnd

                if S.backIndex == S.middleStart:
                    # Only happens if you pop from an empty list - fail nicely
                    S.backIndex -= 1

            # In case of length 0, when everything points at this value
            S.buffer[S.frontIndex & S.bufferMask] = PeakHold.lowest

        S.backIndex += 1
        if S.workingIndex != S.middleStart:
            S.workingIndex -= 1
            S.buffer[S.workingIndex & S.bufferMask] = S.workingMax = max(
                [S.workingMax, S.buffer[S.workingIndex & S.bufferMask]]
            )

    def read(S):
        backMax: float = S.buffer[S.backIndex & S.bufferMask]
        return max([backMax, S.middleMax, S.frontMax])

    def process(S, v: float):
        S.push(v)
        S.pop()
        return S.read()


def peakHold2(sig, holdTime, neutral=0):
    out = np.zeros_like(sig)
    ph = PeakHold(holdTime)
    for i in range(len(sig)):
        out[i] = ph.process(sig[i])
    return out


def peakHoldTarget(sig, holdTime, neutral=0):
    out = np.zeros_like(sig)
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


if __name__ == "__main__":
    # samplerate = 48000
    # rng = np.random.default_rng(3)
    # sig = pulseNoise(rng, 1 / 8, 256)
    sig, fs = soundfile.read("snd/input.wav")

    holdtime = 31
    delay = holdtime - 1
    th = 0.05
    triangleFir = makeTriangleFir(delay)

    holdedTarget = peakHoldTarget(np.abs(sig), holdtime)
    holdedSmithPy = peakHold2(np.abs(sig), holdtime)
    holdedSmithCpp, fs = soundfile.read("snd/output_peak.wav")
    delayed = np.hstack((np.zeros(delay), sig[:-delay]))

    smoothedTarget = signal.lfilter(triangleFir, 1, holdedTarget)
    smoothedSmithPy = signal.lfilter(triangleFir, 1, holdedSmithPy)
    smoothedSmithCpp = signal.lfilter(triangleFir, 1, holdedSmithCpp)

    plt.plot(delayed, alpha=0.5, color="purple", lw=1)
    # plt.plot(smoothedTarget, label="A")
    # plt.plot(smoothedSmithPy, label="B")
    plt.plot(sig, alpha=0.5, color="black", lw=1)
    plt.plot(holdedTarget, label="Target", color="blue", alpha=0.5, lw=4)
    plt.plot(holdedSmithPy, label="Smith Py.", color="orange", alpha=0.75, lw=4)
    plt.plot(holdedSmithCpp, label="Smith C++", color="red", alpha=0.75, lw=1)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()
