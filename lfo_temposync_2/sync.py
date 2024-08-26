import numpy as np
import matplotlib.pyplot as plt


def generatePhase(normalizedFrequency, targetPhase: float):
    sig = np.zeros_like(normalizedFrequency)
    phase = targetPhase - np.floor(targetPhase)
    for i in range(len(sig)):
        sig[i] = phase
        if i >= 1 and normalizedFrequency[i] != normalizedFrequency[i - 1]:
            phase += 0.5
        phase += normalizedFrequency[i]
        phase -= np.floor(phase)
    return sig


def triangle(phase):
    phi = phase + 0.25
    phi -= np.floor(phi)
    return 1 - 4 * np.abs(phi - 0.5)


def saw(phase, pw):
    phi = phase - (1 - pw)
    return np.where(phi < 0, -phi / (1 - pw), phi / pw)


def square(phase):
    return np.where(phase < 0.5, 1, -1)


def crudeSin(phase):
    if phase < 0.25:
        x = 1 - 4 * phase
        return 1 - x * x
    elif phase < 0.5:
        x = 4 * (phase - 0.25)
        return 1 - x * x
    elif phase < 0.75:
        x = 1 - 4 * (phase - 0.5)
        return x * x - 1
    else:
        x = 4 * (phase - 0.75)
        return x * x - 1


def testOsc():
    phase = np.linspace(0, 1, 1024)
    sig = np.zeros_like(phase)
    for i in range(len(sig)):
        sig[i] = triangle(phase[i])
    plt.plot(sig)
    plt.grid()
    plt.show()


def syncKuramoto(targetPhase, initialFrequency, initialPhase, syncRate=0.01):
    """
    Kuramoto model: https://en.wikipedia.org/wiki/Kuramoto_model

    It also works with triangle wave instead of sine.
    It doesn't work when target frequency is too high.
    Phase slightly delays.
    """
    out = np.zeros_like(targetPhase)
    phase = initialPhase - np.floor(initialPhase)
    freq = initialFrequency
    for i in range(len(targetPhase)):
        phase += freq + syncRate * np.sin(2 * np.pi * (targetPhase[i] - phase))
        # phase += freq + syncRate * triangle(targetPhase[i] - phase)
        phase -= np.floor(phase)
        out[i] = phase
    return out


def syncKuramoto2(targetPhase, initialFrequency, initialPhase, syncRate=0.01):
    out = np.zeros_like(targetPhase)
    phase = initialPhase - np.floor(initialPhase)
    freq = initialFrequency
    for i in range(len(targetPhase)):
        out[i] = phase

        p1 = phase
        phase += freq + syncRate * np.sin(2 * np.pi * (targetPhase[i] - phase))
        phase -= np.floor(phase)

        estimatedFreq = phase - p1
        estimatedFreq -= np.floor(estimatedFreq)
        freq += (1 / 2) * syncRate * (estimatedFreq - freq)
    return out


def syncKuramoto3(targetPhase, initialFrequency, initialPhase, syncRate=0.01, nStage=4):
    freqSyncRate = (1 / 64) * syncRate

    out = np.zeros_like(targetPhase)
    phase = np.full(nStage, initialPhase - np.floor(initialPhase))
    freq = initialFrequency
    resultFreq = np.zeros_like(out)

    def increment(target, phase, freq):
        phase += freq + syncRate * np.sin(2 * np.pi * (target - phase))
        phase -= np.floor(phase)
        return phase

    for i in range(len(targetPhase)):
        out[i] = phase[-1]
        resultFreq[i] = freq

        phase[0] = increment(targetPhase[i], phase[0], freq)
        for j in range(1, nStage - 1):
            phase[j] = increment(phase[j - 1], phase[j], freq)

        p1 = phase[-1]
        phase[-1] = increment(phase[-2], phase[-1], freq)

        estimatedFreq = phase[-1] - p1 if phase[-1] >= p1 else phase[-1] - p1 + 1
        freq += freqSyncRate * (estimatedFreq - freq)
    return (out, resultFreq)


def syncFilter(targetPhase, targetFreq, initialFreq, initialPhase, syncRate=0.01):
    out = np.zeros_like(targetPhase)
    phase = initialPhase - np.floor(initialPhase)
    delta = initialFreq
    for i in range(len(targetPhase)):
        phase += delta
        phase -= np.floor(phase)

        d1 = targetPhase[i] - phase
        if d1 < 0:
            d2 = d1 + 1
            phase += syncRate * (d2 if d2 < -d1 else d1)
        else:
            d2 = d1 - 1
            phase += syncRate * (d2 if -d2 < d1 else d1)

        delta += syncRate * (targetFreq[i] - delta)

        out[i] = phase
    return out


def syncFilterAbsBad(targetPhase, targetFreq, initialFreq, initialPhase, syncRate=0.01):
    """This one doesn't work well."""
    dtype = np.float64
    out = np.zeros_like(targetPhase, dtype=dtype)
    phase = (initialPhase - np.floor(initialPhase)).astype(dtype)
    delta = np.array(initialFreq).astype(dtype)
    small = 1 / 2**10
    for i in range(len(targetPhase)):
        phase += delta
        phase -= np.floor(phase)

        d1 = targetPhase[i] - phase
        if d1 < 0:
            d2 = d1 + 1
            phase += syncRate * abs(d2 if d2 < -d1 else d1)
        else:
            d2 = d1 - 1
            phase += syncRate * abs(d2 if -d2 < d1 else d1)

        delta += syncRate * (targetFreq[i] - delta)

        out[i] = phase
    return out


def syncFilterAbs(targetPhase, targetFreq, initialFreq, initialPhase, syncRate=0.01):
    dtype = np.float64
    out = np.zeros_like(targetPhase, dtype=dtype)
    phase = (initialPhase - np.floor(initialPhase)).astype(dtype)
    delta = np.array(initialFreq).astype(dtype)
    small = 1 / 2**10
    for i in range(len(targetPhase)):
        phase += delta
        phase -= np.floor(phase)

        d1 = targetPhase[i] - phase
        diff = (
            (d1 + 1 if d1 + 1 < -d1 else d1)
            if d1 < 0
            else (d1 - 1 if -(d1 - 1) < d1 else d1)
        )
        absed = abs(diff)
        phase += syncRate * (absed if absed >= small else diff)

        delta += syncRate * (targetFreq[i] - delta)

        out[i] = phase
    return out


def plotSync():
    sampleRate = 48000
    targetFreq = 10 / sampleRate
    targetPhase = 0.5
    initialFreq = 50 / sampleRate
    initialPhase = 0.0
    length = 2**15
    syncRate = 0.01

    tgtFreq = np.hstack(
        [
            np.full(length, targetFreq),
            # np.full(length, targetFreq * 4),
            # np.full(length, targetFreq * 20),
        ]
    )
    tgtPhase = generatePhase(tgtFreq, targetPhase)

    # outPhase = syncKuramoto(tgtPhase, initialFreq, initialPhase, syncRate)
    # outPhase = syncKuramoto2(tgtPhase, initialFreq, initialPhase, syncRate)
    # outPhase = syncFilter(tgtPhase, tgtFreq, initialFreq, initialPhase, syncRate)
    # outPhase = syncFilterAbsBad(tgtPhase, tgtFreq, initialFreq, initialPhase, syncRate)
    outPhase = syncFilterAbs(tgtPhase, tgtFreq, initialFreq, initialPhase, syncRate)

    fig = plt.figure()
    fig.set_size_inches((8, 3))
    plt.title(f"Forward Only EMA Filter Synchronization")
    plt.plot(tgtPhase, color="blue", alpha=0.5, lw=1, label="Target")
    plt.plot(outPhase, color="red", alpha=0.5, lw=1, label="Follower")
    # plt.xlabel("Time [sample]")
    plt.ylabel("Phase [rad/2π]")
    plt.legend(loc=10, bbox_to_anchor=(0.66, 0.5))
    plt.grid()
    plt.tight_layout()
    plt.show()


def testSync():
    sampleRate = 48000
    targetFreq = 10 / sampleRate
    targetPhase = 0.5
    initialFreq = 50 / sampleRate
    initialPhase = 0.0
    length = 4096
    syncRate = 0.01

    tgtFreq = np.hstack(
        [
            np.full(length, targetFreq),
            np.full(length, targetFreq * 4),
            np.full(length, targetFreq * 20),
        ]
    )
    tgtPhase = generatePhase(tgtFreq, targetPhase)

    # outPhase = syncKuramoto(tgtPhase, initialFreq, initialPhase, syncRate)
    # outPhase = syncKuramoto2(tgtPhase, initialFreq, initialPhase, syncRate)
    # outPhase = syncFilter(tgtPhase, tgtFreq, initialFreq, initialPhase, syncRate)
    outPhase = syncFilterAbs(tgtPhase, tgtFreq, initialFreq, initialPhase, syncRate)

    phaseError = calcPhaseError(tgtPhase, outPhase)
    freqError = calcFrequencyError(tgtPhase, outPhase)

    fig, ax = plt.subplots(3)
    fig.set_size_inches((8, 8))
    fig.suptitle(f"Forward Only EMA Filter Synchronization on Varying Target Frequency")
    ax[0].plot(tgtPhase, color="blue", alpha=0.5, lw=1, label="Target")
    ax[0].plot(outPhase, color="red", alpha=0.5, lw=1, label="Follower")
    # ax[0].set_xlabel("Time [sample]")
    ax[0].set_ylabel("Phase [rad/2π]")
    ax[0].legend()

    ax[1].plot(phaseError, color="red", alpha=0.5, lw=1)
    # ax[1].set_xlabel("Time [sample]")
    ax[1].set_ylabel("Phase Error [rad/2π]")

    ax[2].plot(freqError, color="red", alpha=0.5, lw=1)
    ax[2].set_xlabel("Time [sample]")
    ax[2].set_ylabel("Frequency Error [rad/(2π * sample)]")
    ax[2].set_ylim((-1e-2, 1e-2))

    for axis in ax:
        axis.grid()
    fig.tight_layout()
    plt.show()


def testSyncWithNoise():
    sampleRate = 48000
    targetFreq = 10 / sampleRate
    targetPhase = 0.5
    initialFreq = 200 / sampleRate
    initialPhase = 0.0
    length = 2**17
    syncRate = 0.01
    noiseGain = 1 / 2
    rng = np.random.default_rng(0)

    tgtFreq = np.hstack(
        [
            np.full(length, targetFreq),
            # np.full(length, targetFreq * 4),
            # np.full(length, targetFreq * 20),
        ]
    )
    tgtPhase = generatePhase(tgtFreq, targetPhase)

    tgtNoise = tgtPhase + rng.normal(0, noiseGain / 3, len(tgtPhase))
    tgtNoise -= np.floor(tgtNoise)

    outPhase, resultFreq = syncKuramoto3(
        tgtNoise, initialFreq, initialPhase, syncRate, 40
    )

    phaseError = calcPhaseError(tgtPhase, outPhase)
    freqError = calcFrequencyError(tgtPhase, outPhase)

    fig, ax = plt.subplots(4)
    fig.set_size_inches((8, 9))
    fig.suptitle(
        f"Nested Kuramoto Model Synchronization on Noisy Signal (S/N={int(1/noiseGain)}/1)"
    )
    ax[0].plot(tgtNoise, color="black", alpha=0.25, lw=1, label="Target + Noise")
    ax[0].plot(tgtPhase, color="blue", alpha=0.5, lw=1, label="Target")
    ax[0].plot(outPhase, color="red", alpha=0.5, lw=1, label="Follower")
    # ax[0].set_xlabel("Time [sample]")
    ax[0].set_ylabel("Phase [rad/2π]")
    # ax[0].set_ylim((-0.5, 1.5))
    ax[0].legend()

    ax[1].plot(phaseError, color="red", alpha=0.5, lw=1)
    # ax[1].set_xlabel("Time [sample]")
    ax[1].set_ylabel("Phase Error")

    ax[2].plot(freqError, color="red", alpha=0.5, lw=1)
    ax[2].set_ylabel("Freq. Error")
    # ax[2].set_xlabel("Time [sample]")
    ax[2].set_ylim((-1e-2, 1e-2))

    ax[3].plot(tgtFreq, color="blue", alpha=0.5, lw=1, label="Target")
    ax[3].plot(resultFreq, color="red", alpha=0.5, lw=1, label="Follower")
    ax[3].set_ylabel("Follower Freq. [rad/2π]")
    ax[3].set_xlabel("Time [sample]")
    ax[3].legend()

    for axis in ax:
        axis.grid()
    fig.tight_layout()
    plt.show()


def calcPhaseError(target, follower):
    diff = follower - target
    diff[diff < -0.5] += 1
    diff[diff > 0.5] -= 1
    return diff

    # # debug
    # plt.plot(diff, color="black", alpha=0.5)
    # plt.plot(follower, color="red", alpha=0.5)
    # plt.plot(target, color="blue", alpha=0.5)
    # plt.grid()
    # plt.show()


def calcFrequencyError(target, follower):
    return np.diff(np.unwrap(follower - target, period=1))


def plotErrors():
    sampleRate = 48000
    targetFreq = np.geomspace(1 / 48000, 1 / 2**2, 16)
    targetPhase = 0
    initialFreq = 1000 / sampleRate
    initialPhase = 0.5
    length = 2048
    syncRate = 0.1

    cmap = plt.get_cmap("magma")
    fig, ax = plt.subplots(2)
    for index, rFrq in enumerate(targetFreq):
        tgtFreq = np.full(length, rFrq)
        tgtPhase = generatePhase(tgtFreq, targetPhase)

        outPhase = syncKuramoto(tgtPhase, initialFreq, initialPhase, syncRate)
        # outPhase = syncFilter(tgtPhase, tgtFreq, initialFreq, initialPhase, syncRate)

        phaseError = calcPhaseError(tgtPhase, outPhase)
        freqError = calcFrequencyError(tgtPhase, outPhase)
        ax[0].plot(
            phaseError,
            color=cmap(index / len(targetFreq)),
            lw=1,
            alpha=0.5,
        )
        ax[1].plot(
            freqError,
            color=cmap(index / len(targetFreq)),
            lw=1,
            alpha=0.5,
        )
    ax[0].set_ylabel("Phase Error")
    ax[1].set_ylabel("Frequency Error")
    for axis in ax:
        axis.set_xlabel("Time [sample]")
        axis.grid()
        # axis.legend()
    fig.set_size_inches((8, 8))
    fig.tight_layout()
    plt.show()


def plotErrorVsFrequency():
    sampleRate = 48000
    targetFreq = np.geomspace(1 / 2**12, 1 / 2, 256)
    targetPhase = 0
    initialFreq = 0.002
    initialPhase = 0.5
    length = 8192
    syncRate = [0.1, 0.01, 0.001, 0.0001]

    fig, ax = plt.subplots(2)
    fig.set_size_inches(6, 6)
    fig.suptitle(f"Synchronization Error(EMA Filter, ω={initialFreq:.3f})")
    cmap = plt.get_cmap("hot")

    for index, syncRt in enumerate(syncRate):
        phaseError = []
        freqError = []
        for rFrq in targetFreq:
            tgtFreq = np.full(length, rFrq)
            tgtPhase = generatePhase(tgtFreq, targetPhase)

            outPhase = syncKuramoto(tgtPhase, initialFreq, initialPhase, syncRt)
            # outPhase = syncFilter(tgtPhase, tgtFreq, initialFreq, initialPhase, syncRt)

            half = len(tgtPhase) // 2
            d0 = calcPhaseError(tgtPhase, outPhase)
            phaseError.append(np.average(np.abs(d0[half:])))
            d1 = calcFrequencyError(tgtPhase, outPhase)
            freqError.append(np.average(np.abs(d1[half:])))

        freqHz = targetFreq
        color = cmap(index / len(syncRate))
        ax[0].plot(freqHz, phaseError, color=color, alpha=0.75, label=f"K={syncRt}")
        ax[1].plot(freqHz, freqError, color=color, alpha=0.75, label=f"K={syncRt}")

    initHz = initialFreq
    ax[0].axvline(initHz, color="black", lw=1, alpha=0.75, ls="--", label="ω")
    ax[0].set_ylabel("Phase Error Avg. [rad/2π]")
    ax[0].set_xlabel("Target Frequency [rad/2π]")

    ax[1].axvline(initHz, color="black", lw=1, alpha=0.75, ls="--", label="ω")
    ax[1].set_ylabel("Frequency Error Avg. [rad/2π]")
    ax[1].set_xlabel("Target Frequency [rad/2π]")
    # ax[1].set_yscale("log")
    freqErrorRange = 0.01
    freqErrorMargin = freqErrorRange / 10
    ax[1].set_ylim((-freqErrorMargin, freqErrorRange + freqErrorMargin))

    for axis in ax:
        axis.set_xscale("log")
        axis.grid(which="both", color="#f0f0f0")
        axis.legend()
    fig.tight_layout()
    plt.show()


def kuramotoOrderParameterR():
    """
    Plot order parameter `r` in following paper.

    Reference: Daniels, Bryan C. "Synchronization of globally coupled nonlinear oscillators: the rich behavior of the Kuramoto model." Ohio Wesleyan Physics Dept., Essay 7, no. 2 (2005): 20.
    """
    rng = np.random.default_rng()
    frequency = np.pi * 2 ** rng.uniform(-4, 0, 4)

    time = np.linspace(0, 10, 2**16)
    summed = np.zeros_like(time, dtype=np.complex128)
    cmap = plt.get_cmap("viridis")
    for index, frq in enumerate(frequency):
        sig = np.exp(1j * frq * time)
        summed += sig
        color = cmap(index / len(frequency))
        plt.plot(time, sig.imag, color=color, alpha=0.5, lw=1, label=index)
    r = np.abs(summed) / len(frequency)

    plt.plot(time, r, color="black", label="r")
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # testOsc()
    # plotSync()
    # testSync()
    # testSyncWithNoise()
    # plotErrors()
    # plotErrorVsFrequency()
    kuramotoOrderParameterR()
