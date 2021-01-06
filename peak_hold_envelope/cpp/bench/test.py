import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile

from collections import deque

def plotSignal(signals, labels, title, legendLocation=0):
    plt.figure(figsize=(6, 3), tight_layout=True)
    plt.title(title)
    cmap = plt.get_cmap("viridis")
    for idx, (sig, label) in enumerate(zip(signals, labels)):
        plt.plot(sig, lw=1, color=cmap(idx / len(signals)), label=label)
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.legend(loc=legendLocation)
    plt.grid()

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

def plot():
    dataIdeal, _ = soundfile.read("snd/Ideal.wav", always_2d=True)
    outI = dataIdeal.T[0]

    data, samplerate = soundfile.read("snd/input.wav", always_2d=True)
    sig = data.T[0]

    holdI = idealPeakHoldFast(sig, 64)

    plotSignal([sig, holdI, outI], ["input", "Ref.", "C++"], "Ideal Peak Hold")
    plt.show()

def testDiff(path1, path2):
    data1, fs1 = soundfile.read(path1, always_2d=True)
    data2, fs2 = soundfile.read(path2, always_2d=True)
    assert fs1 == fs2
    np.testing.assert_almost_equal(data1, data2)

def test(holdtime=64):
    print("Testing the output of C++ implementation.")

    testDiff("snd/IdealDeque.wav", "snd/Ideal.wav")

    sig, _ = soundfile.read("snd/input.wav", always_2d=False)

    dataI, _ = soundfile.read("snd/Ideal.wav", always_2d=False)
    holdI = idealPeakHoldFast(sig, holdtime)
    np.testing.assert_almost_equal(dataI, holdI)

if __name__ == "__main__":
    test(64)
