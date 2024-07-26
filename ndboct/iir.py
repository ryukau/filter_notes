import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def plotReferenceSlope(startHz, slopeDB, step, nCascade):
    refX = []
    refY = []
    for idx in range(0, step * nCascade, step):
        refX.append(startHz * 2**idx)
        refY.append(idx * slopeDB)
    plt.plot(refX, refY, label=f"Ref. {slopeDB} dB/oct", color="red", ls="--")


def toDecibel(response):
    return 20 * np.log10(np.maximum(np.abs(response), 1e-7))


def onepoleHighShelfEma(cutoffNormalized, gainDB):
    """EMA stands for exponential moving average."""
    y = 1 - np.cos(2 * np.pi * cutoffNormalized)
    k = np.sqrt(y * y + 2 * y) - y
    g = 10 ** (gainDB / 20)

    a0 = 1
    a1 = k - 1
    a2 = 0
    b0 = (1 - g) * k + g * a0
    b1 = (1 - g) * 0 + g * a1
    b2 = 0
    return [b0 / a0, b1 / a0, b2 / a0, 1, a1 / a0, a2 / a0]


def onepoleHighShelfBilinear(cutoffNormalized, gainDB):
    k = 1 / np.tan(np.pi * cutoffNormalized)
    g = 10 ** (gainDB / 20)

    a0 = 1 + k
    a1 = 1 - k
    a2 = 0
    b0 = (1 - g) + g * a0
    b1 = (1 - g) + g * a1
    b2 = 0
    return [b0 / a0, b1 / a0, b2 / a0, 1, a1 / a0, a2 / a0]


def onepoleHighShelfMatched(cutoffNormalized, gainDB):
    fc = 2 * cutoffNormalized
    G = 10 ** (gainDB / 20)

    fm = 0.9
    φm = 1 - np.cos(np.pi * fm)

    def alphabeta(V):
        return 2 / (np.pi * np.pi) * (1 / (fm * fm) + V) - 1 / φm

    α = alphabeta(1 / (G * fc * fc))
    β = alphabeta(G / (fc * fc))

    a1 = -α / (1 + α + np.sqrt(1 + 2 * α))
    b = -β / (1 + β + np.sqrt(1 + 2 * β))
    b0 = (1 + a1) / (1 + b)
    b1 = b * b0

    return [b0, b1, 0, 1, a1, 0]


def testOnepoleHighShelf():
    sampleRate = 48000
    cutoff = 1000 / sampleRate
    gain = 20 * np.log10(0.5)

    sosEma = onepoleHighShelfEma(cutoff, gain)
    sosBlt = onepoleHighShelfBilinear(cutoff, gain)
    sosMtd = onepoleHighShelfMatched(cutoff, gain)

    worN = np.geomspace(1, sampleRate / 2, 2**12)

    cmap = plt.get_cmap("cool")

    f, r = signal.sosfreqz(sosEma, worN=worN, fs=sampleRate)
    plt.plot(f, toDecibel(r), color=cmap(0 / 3), alpha=0.5, label="EMA")

    f, r = signal.sosfreqz(sosBlt, worN=worN, fs=sampleRate)
    plt.plot(f, toDecibel(r), color=cmap(1 / 3), alpha=0.5, label="Bilinear")

    f, r = signal.sosfreqz(sosMtd, worN=worN, fs=sampleRate)
    plt.plot(f, toDecibel(r), color=cmap(2 / 3), alpha=0.5, label="Matched")

    plt.xscale("log")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Gain [dB]")
    plt.ylim((-12, 2))
    plt.grid()
    plt.legend()
    plt.show()


def doubleOrder(sosX, sosY):
    """Only works for 1-pole filter."""
    sos = np.ones(6)

    x_b0 = sosX[0]
    x_b1 = sosX[1]
    y_b0 = sosY[0]
    y_b1 = sosY[1]
    sos[0] = x_b0 * y_b0
    sos[1] = x_b0 * y_b1 + y_b0 * x_b1
    sos[2] = x_b1 * y_b1

    x_a1 = sosX[4]
    y_a1 = sosY[4]
    sos[4] = x_a1 + y_a1
    sos[5] = x_a1 * y_a1

    return sos


def testOnePoleSlope(onepoleFuncs):
    """
    onepoleFuncs = [
        {func: somefunction, label: "label"},
        # ...
    ]
    """
    sampleRate = 48000
    startHz = 20
    slope = -1
    nCascade = 10

    cmap = plt.get_cmap("cool")
    for index, flt in enumerate(onepoleFuncs):
        sos = []
        cutoffHz = startHz
        for _ in range(nCascade):
            sos.append(flt["func"](cutoffHz / sampleRate, slope))

            # # Below 2 lines are an experimental implementation.
            # sos.append(flt["func"](cutoffHz / sampleRate, slope / 2))
            # sos[-1] = doubleOrder(sos[-1], sos[-1])

            cutoffHz *= 2
        freq, response = signal.sosfreqz(
            sos, worN=np.geomspace(1, sampleRate / 2, 2**16), fs=sampleRate
        )
        plt.plot(
            freq,
            toDecibel(response),
            label=flt["label"],
            color=cmap(index / len(onepoleFuncs)),
            lw=2,
        )

    plotReferenceSlope(startHz, slope, 2, nCascade // 2)

    plt.title(f"Gain Response of Slope Filter ({slope} dB/oct)")
    plt.xscale("log")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Gain [dB]")
    plt.grid(which="both", color="#f4f4f4")
    plt.legend()
    plt.tight_layout()
    plt.show()


def biquadHighShelf(cutoffNormalized, Q, gainDB):
    ω0 = 2 * np.pi * cutoffNormalized
    cs = np.cos(ω0)
    sn = np.sin(ω0)
    A = np.power(10, gainDB / 40)
    α = 0.5 * sn * np.sqrt((A + 1 / A) * (1 / Q - 1) + 2)

    B = 2 * np.sqrt(A) * α

    a0 = (A + 1) - (A - 1) * cs + B
    a1 = 2 * ((A - 1) - (A + 1) * cs)
    a2 = (A + 1) - (A - 1) * cs - B
    b0 = A * ((A + 1) + (A - 1) * cs + B)
    b1 = -2 * A * ((A - 1) + (A + 1) * cs)
    b2 = A * ((A + 1) + (A - 1) * cs - B)

    # print((-a1 + np.sqrt(a1 * a1 - 4 * a0 * a2 + 0j)) / (2 * a0))

    return [b0 / a0, b1 / a0, b2 / a0, 1, a1 / a0, a2 / a0]


def getSlopeSos(sampleRate, startHz, Q, slopeDB, step, nCascade):
    sos = []
    cutoffHz = startHz * np.sqrt(2 * step)
    for _ in range(0, step * nCascade, step):
        sos.append(biquadHighShelf(cutoffHz / sampleRate, Q, step * slopeDB))
        cutoffHz *= 2 * step
    return sos


def plotSos(name, sos, sampleRate, color):
    freq, response = signal.sosfreqz(
        sos, worN=np.geomspace(1, sampleRate / 2, 2**16), fs=sampleRate
    )
    plt.plot(freq, toDecibel(response), label=name, color=color, lw=2)
    # plt.plot(freq, toPhase(response), label=name, color=color, lw=4)
    # plt.plot(freq[:-1], toGroupDelay(response), label=name, color=color, lw=4)
    # plt.plot(np.diff(toDecibel(response)), label=name, color=color)


def testIIRSlope():
    sampleRate = 48000
    startHz = 20
    slope = -1
    nCascade = 10

    cmap = plt.get_cmap("hot")

    maxStep = 4
    for step in range(1, maxStep + 1):
        Q = np.sqrt(2) / 2
        sos = getSlopeSos(sampleRate, startHz, Q, slope, step, nCascade // step)
        plotSos(f"k={step}", sos, sampleRate, cmap((step - 1) / maxStep))

    plotReferenceSlope(startHz * 1.25, slope, 2, nCascade // 2)

    plt.title(f"Gain Response of Slope Filter ({slope} dB/oct)")
    plt.xscale("log")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Gain [dB]")
    plt.xscale("log")
    plt.grid(which="both", color="#f4f4f4")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # testOnepoleHighShelf()
    # testOnePoleSlope(
    #     [
    #         {"func": onepoleHighShelfEma, "label": "EMA"},
    #         {"func": onepoleHighShelfBilinear, "label": "Bilinear"},
    #         {"func": onepoleHighShelfMatched, "label": "Matched"},
    #     ]
    # )
    testIIRSlope()
