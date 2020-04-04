import numpy as np
import scipy.signal as signal
import scipy.optimize as optimize
import matplotlib.pyplot as pyplot

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

def plot(c=-0.01, k=0.96):
    cmap = pyplot.get_cmap("viridis")
    b, a = transferFunction(-0.2, 0.9)
    w, h = signal.freqz(b, a, 2**16)
    pyplot.plot(
        w,
        20 * np.log10(abs(h)),
        color=cmap(0.5),
        label=f"c={c:.3f},k={k:.3f}",
    )
    pyplot.ylabel("Amplitude [dB]", color="b")
    pyplot.xlabel("Frequency [rad/sample]")
    pyplot.xscale("log")
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

def cutoffPlot():
    cmap = pyplot.get_cmap("viridis")
    nPlot = 12
    for idx, c in enumerate(np.geomspace(1e-5, 1, nPlot)):
        b, a = transferFunction(c, 0.0)
        ω, h = signal.freqz(b, a, 2**16)
        gain = 20 * np.log10(abs(h))
        pyplot.plot(ω, gain, color=cmap(idx / nPlot), label=f"c={c:.3f}")
        index = np.argmax(gain <= -3)
        if index < 2:
            continue
        prev = index - 1
        delta_range = gain[index] - gain[prev]
        ratio = (-3 - gain[prev]) / delta_range
        x = ω[prev] + ratio * (ω[index] - ω[prev])
        y = gain[prev] + ratio * delta_range
        pyplot.scatter(x, y, color=cmap(idx / nPlot))
    pyplot.title("Amplitude Response, k=0")
    pyplot.ylabel("Amplitude [dB]", color="b")
    pyplot.xlabel("Frequency [rad/sample]")
    pyplot.xscale("log")
    pyplot.legend()
    pyplot.grid()
    pyplot.show()

def highpassPlot():
    cmap = pyplot.get_cmap("viridis")
    nPlot = 32
    for idx, α in enumerate(1 + 1e-5 - np.geomspace(1e-5, 0.5, nPlot)):
        b, a = transferFunction(1, 0, α)
        ω, h = signal.freqz(b, a, 2**16)
        gain = 20 * np.log10(abs(h))
        pyplot.plot(ω, gain, color=cmap(idx / nPlot), label=f"α={α:.5f}")

        index = np.argmax(gain >= -3)
        if index < 2:
            continue
        prev = index - 1
        delta_range = gain[index] - gain[prev]
        ratio = (-3 - gain[prev]) / delta_range
        x = ω[prev] + ratio * (ω[index] - ω[prev])
        y = gain[prev] + ratio * delta_range
        pyplot.scatter(x, y, color=cmap(idx / nPlot))
    pyplot.ylabel("Amplitude [dB]", color="b")
    pyplot.xlabel("Frequency [rad/sample]")
    pyplot.xscale("log")
    pyplot.grid()
    pyplot.legend(loc=4)
    pyplot.show()

def dampingCutoffPlot():
    cmapReal = pyplot.get_cmap("viridis")
    cmapEstimate = pyplot.get_cmap("plasma")
    xC = np.linspace(1e-5, 1, 128)
    data = []
    for idx, c in enumerate(xC):
        b, a = transferFunction(c, 0.0)
        ω, h = signal.freqz(b, a, 2**16)
        gain = 20 * np.log10(abs(h))
        index = np.argmax(gain <= -3)
        if index < 2:
            continue
        prev = index - 1
        delta_range = gain[index] - gain[prev]
        ratio = (-3 - gain[prev]) / delta_range
        cutoff = ω[prev] + ratio * (ω[index] - ω[prev])
        data.append((c, cutoff))
        pyplot.scatter(cutoff, c, color=cmapReal(idx / len(xC)))

    xC_fit, ω_c_fit = zip(*data)
    ω_c_fit = np.array(ω_c_fit) / 2 / np.pi
    polyCoef = np.polyfit(ω_c_fit, xC_fit, 6)
    print(ω_c_fit)
    print([p for p in polyCoef])
    polyCoef[-1] = 0
    poly = np.poly1d(polyCoef)
    curveX = np.linspace(0, 0.5, 1024)
    curveY = poly(curveX)
    pyplot.plot(curveX * 2 * np.pi, curveY, label="polyfit")

    pyplot.ylabel("Damping Factor")
    pyplot.xlabel("Cutoff Frequency [rad/sample]")
    # pyplot.xscale("log")
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

def highpassCutoffPlot():
    cmap = pyplot.get_cmap("viridis")
    x_α = 1 + 1e-5 - np.linspace(1e-5, 1, 128)
    data = []
    for idx, α in enumerate(x_α):
        b, a = transferFunction(1, 0, α)
        ω, h = signal.freqz(b, a, 2**16)
        gain = 20 * np.log10(abs(h))
        index = np.argmax(gain >= -12)
        if index < 2:
            continue
        prev = index - 1
        delta_range = gain[index] - gain[prev]
        ratio = (-3 - gain[prev]) / delta_range
        cutoff = ω[prev] + ratio * (ω[index] - ω[prev])
        data.append((α, cutoff))
        pyplot.scatter(cutoff, α, color=cmap(idx / len(x_α)))

    x_α_fit, ω_c_fit = zip(*data)

    def curve_func(x, p0, p1):
        return (1 - p1) + p1 * np.exp(p0 * x)

    popt, _ = optimize.curve_fit(
        curve_func, np.array(ω_c_fit) / np.pi, x_α_fit, bounds=([-np.inf, 0], [0, 1]))
    print([p for p in popt])
    plotX = np.linspace(0, 1, 1024)
    curve = curve_func(plotX, *popt)
    pyplot.plot(plotX * np.pi, curve, label="expfit")

    x_α_fit, ω_c_fit = zip(*data)
    ω_c_fit = np.array(ω_c_fit) / 2 / np.pi
    polyCoef = np.polyfit(ω_c_fit, x_α_fit, 8)
    print([p for p in polyCoef])
    poly = np.poly1d(polyCoef)
    curveX = np.linspace(0, 0.5, 1024)
    curveY = poly(curveX)
    pyplot.plot(curveX * 2 * np.pi, curveY, label="polyfit")

    pyplot.title("α-ω_c Plot")
    pyplot.ylabel("α")
    pyplot.xlabel("Cutoff Frequency [rad/sample]")
    # pyplot.xscale("log")
    pyplot.xlim([-0.1, np.pi + 0.1])
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

def resonanceFrequencyPlot():
    cmap = pyplot.get_cmap("viridis")
    nC = 8
    nK = 128
    for idx, c in enumerate(np.linspace(1e-3, 1, nC)):
        dataX = []
        dataY = []
        for k in (1 - np.geomspace(1e-5, 1, nK)):
            b, a = transferFunction(c, k)
            ω, h = signal.freqz(b, a, 2**16)
            h[np.isfinite(h) == False] = 0
            gain = 20 * np.log10(abs(h))
            index = np.argmax(gain)
            dataX.append(ω[index])
            dataY.append(gain[index])
        pyplot.plot(dataX, dataY, color=cmap(idx / nC), label=f"c={c:.3f}")
    pyplot.ylabel("Peak Gain [dB]")
    pyplot.xlabel("Cutoff Frequency [rad/sample]")
    # pyplot.xscale("log")
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

def resonancePeakPlot():
    cmap = pyplot.get_cmap("viridis")
    nC = 10
    nK = 128
    kk = 1 - np.geomspace(1e-5, 1, nK)
    for idx, c in enumerate(np.linspace(0.1, 1, nC)):
        dataY = []
        for k in kk:
            b, a = transferFunction(c, k, 1)
            ω, h = signal.freqz(b, a, 2**16)
            h[np.isfinite(h) == False] = 0
            peak = np.max(np.abs(h))
            if not np.isfinite(peak):
                print(idx, c, k, peak)
            dataY.append(peak)
        pyplot.plot(1 - kk, dataY, color=cmap(idx / nC), label=f"c={c:.3f}")
    pyplot.ylabel("Peak Amplitude")
    pyplot.xlabel("1 - k")
    pyplot.xscale("log")
    pyplot.yscale("log")
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

plot(0.1411257617811452, 0.998801682380008)
cutoffPlot()
highpassPlot()
highpassCutoffPlot()
dampingCutoffPlot()
resonanceFrequencyPlot()
resonancePeakPlot()
