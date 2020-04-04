import numpy as np
import scipy.signal as signal
import scipy.optimize as optimize
import matplotlib.pyplot as pyplot
import json

def transferFunction(c, k, α=1):
    return (
        [1, -(k + 1), k, 0],
        [
            -(k - 1) / (c * α),
            (k * k - 1) / (c * α) - (k - 1) / α + (k - 1) / c,
            k * (k - 1) / (c * α) - (k * k - 1) / c + k - 1,
            k * (k - 1) / c,
        ],
    )

def findUniformResonance():
    maxIteration = 256

    data = []

    targetPeak = np.geomspace(1, 1e5, 16)
    nC = 128
    xC = np.geomspace(1e-4, 1, nC)

    for target in targetPeak:
        peakData = []
        for c in xC:
            jdx = 0
            diff = 1
            k = 0.5
            delta = 0.5
            peak = 0
            while diff > 1e-5 and jdx < maxIteration:
                b, a = transferFunction(c, k)
                ω, h = signal.freqz(b, a, 2**16)
                h[np.isfinite(h) == False] = 0
                peak = np.max(np.abs(h))
                diff = abs(target - peak)
                delta *= 0.5
                if peak < target:
                    k += delta
                else:
                    k -= delta
                jdx += 1
            print(f"{target}, {c}, {k}, {peak}")
            peakData.append({"c": c, "k": k, "peak": peak})
        data.append({
            "targetPeak": target,
            "c": [d["c"] for d in peakData],
            "k": [d["k"] for d in peakData],
            "peak": [d["peak"] for d in peakData],
        })
    with open("data.json", "w") as fi:
        json.dump(data, fi, indent=2)

def plotResonance(data):
    cmapV = pyplot.get_cmap("viridis")
    cmapP = pyplot.get_cmap("plasma")
    for idx, dat in enumerate(data):
        if dat["targetPeak"] == 1.0:
            continue
        damping = np.array(dat["c"])
        spring = np.array(dat["k"])
        peak = np.array(dat["peak"])
        pyplot.plot(
            damping,
            1 - spring,
            lw=1,
            color=cmapV(idx / len(data)),
            label=f"Peak={dat['targetPeak']:g}")
        kMin = np.min(spring)
        kMax = np.max(spring)
        y = kMax - (kMax - kMin) * np.arccos(1 - np.array(damping)) / (np.pi / 2)
        pyplot.plot(damping, 1 - y, lw=1, ls=":", color=cmapP(idx / len(data)))
    pyplot.title("c-k Plot (dotted line is approximation)")
    pyplot.xlabel("c")
    pyplot.ylabel("1 - k")
    pyplot.yscale("log")
    pyplot.legend(ncol=2)
    pyplot.grid()
    pyplot.tight_layout()
    pyplot.show()

def plotSpringMinMax(data):
    kMin = []
    kMax = []
    target = []
    for idx, dat in enumerate(data):
        damping = dat["c"]
        spring = dat["k"]
        peak = dat["peak"]
        target.append(dat["targetPeak"])
        kMin.append(np.min(spring))
        kMax.append(np.max(spring))
    #
    # kMin plot.
    plotX = np.linspace(0, 1, 1024)
    alpha = np.log(1 - kMin[-1])
    y = 1 - np.exp(alpha * plotX)
    print(alpha)
    pyplot.plot(plotX, y, zorder=3, alpha=0.75, color="#00ff00", label="approx.")
    fitX = np.linspace(0, 1, len(kMin))
    pyplot.plot(fitX, kMin, lw=4, zorder=2, color="black", label="kMin")
    pyplot.title("k_Min-Resonance Plot")
    pyplot.xlabel("Resonance")
    pyplot.ylabel("k")
    pyplot.grid()
    pyplot.legend()
    pyplot.show()
    #
    # KMax plot.
    plotX = np.linspace(0, 1, 1024)
    offset = 1 - kMin[-1]
    alpha = np.log(offset)
    y = kMax[-1] - 0.01 * (np.exp(alpha * plotX) - offset)
    pyplot.plot(plotX, y, zorder=3, alpha=0.75, color="#00ff00", label="approx.")
    print(kMax[-1], alpha, offset)
    pyplot.plot(fitX, kMax, lw=4, zorder=2, color="black", label="kMax")
    pyplot.title("k_Max-Resonance Plot")
    pyplot.xlabel("Resonance")
    pyplot.ylabel("k")
    pyplot.grid(zorder=1)
    pyplot.legend()
    pyplot.tight_layout()
    pyplot.show()

with open("data.json", "r") as fi:
    data = json.load(fi)

findUniformResonance()
plotResonance(data)
plotSpringMinMax(data)
