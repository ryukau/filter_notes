import numpy as np
import scipy.signal as signal
import scipy.optimize as optimize
import matplotlib.pyplot as pyplot

def transferFunctionU1(k1, k2):
    return (
        [0, 1, -1, 0, 0],
        [
            1 / k2,
            (k1 - 3) / k2 + 2,
            (3 - k1) / k2 + k1 - 4,
            -1 / k2 + 2,
        ],
    )

def transferFunctionU2(k1, k2):
    return (
        [
            1,
            k2 + k1 - 3,
            -2 * k2 - k1 + 3,
            k2 - 1,
            0,
        ],
        [
            1 / k2,
            (k1 - 4) / k2,
            (-2 * k1 + 6) / k2 - k1,
            (k1 - 4) / k2 + k1,
            1 / k2,
        ],
    )

def plot(transferFunction, k1, k2):
    cmap = pyplot.get_cmap("viridis")
    b, a = transferFunction(k1, k2)
    w, h = signal.freqz(b, a, 2**16)
    pyplot.plot(w, 20 * np.log10(abs(h)), color=cmap(0.5))
    pyplot.title(f"k1={k1:.3}, k2={k2:.3f}")
    pyplot.ylabel("Amplitude [dB]", color="b")
    pyplot.xlabel("Frequency [rad/sample]")
    pyplot.xscale("log")
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

def cutoffPlot(transferFunction):
    cmap = pyplot.get_cmap("viridis")
    nPlot = 16
    k1 = np.pi
    for idx, k2 in enumerate(np.geomspace(1e-5, 1, nPlot)):
        b, a = transferFunction(k1, k2)
        ω, h = signal.freqz(b, a, 2**16)
        gain = 20 * np.log10(abs(h))
        gain[np.isfinite(gain) == False] = 0
        pyplot.plot(ω, gain, color=cmap(idx / nPlot), label=f"k2={k2:.3f}")
        index = np.argmax(gain[2:] <= -3)
        if index < 1:
            continue
        index += 2
        prev = index - 1
        delta_range = gain[index] - gain[prev]
        ratio = (-3 - gain[prev]) / delta_range
        cutoff = ω[prev] + ratio * (ω[index] - ω[prev])
        y = gain[prev] + ratio * delta_range
        pyplot.scatter(cutoff, y, color=cmap(idx / nPlot))
    pyplot.title(f"k1={k1}")
    pyplot.ylabel("Gain [dB]")
    pyplot.xlabel("Frequency [rad/sample]")
    pyplot.xscale("log")
    pyplot.ylim((-40, 5))
    pyplot.grid()
    pyplot.legend(ncol=3)
    pyplot.show()

def cutoffK2Plot(transferFunction):
    pyplot.figure(figsize=(6, 3))
    cmap = pyplot.get_cmap("viridis")
    nPlot = 128
    k1 = np.pi

    data = []
    for idx, k2 in enumerate(np.geomspace(1e-5, 1, nPlot)):
        b, a = transferFunction(k1, k2)
        ω, h = signal.freqz(b, a, 2**16)
        gain = 20 * np.log10(abs(h))
        gain[np.isfinite(gain) == False] = 0
        index = np.argmax(gain[2:] <= -3)
        if index < 1:
            continue
        index += 2
        prev = index - 1
        delta_range = gain[index] - gain[prev]
        ratio = (-3 - gain[prev]) / delta_range
        cutoff = ω[prev] + ratio * (ω[index] - ω[prev])
        y = gain[prev] + ratio * delta_range
        pyplot.scatter(cutoff, k2, color=cmap(idx / nPlot))

        data.append((cutoff, k2))

    ω_c_fit, k2_fit = zip(*data)

    def curve_func(x, p0, p1):
        return p0 * x + p1 * x * x

    popt, _ = optimize.curve_fit(curve_func, np.array(ω_c_fit) / 2 / np.pi, k2_fit)
    print([p for p in popt])
    print(ω_c_fit[-1], k2_fit[-1])
    plotX = np.linspace(0, 0.125, 2**16)
    curve = curve_func(plotX, *popt)
    curve -= curve[0]
    pyplot.plot(plotX * 2 * np.pi, curve, label="fit")

    pyplot.title(f"k1={k1}")
    pyplot.ylabel("k2")
    pyplot.xlabel("Frequency [rad/sample]")
    # pyplot.xscale("log")
    pyplot.grid()
    pyplot.legend()
    pyplot.tight_layout()
    pyplot.show()

# plot(transferFunctionU2, np.pi, 0.6)
# cutoffPlot(transferFunctionU2)
cutoffK2Plot(transferFunctionU2)
"""
# k2-cutoff approximation for u2
ω_c <= 0.7250965101345954
cutoffHz <~

x = cutoffHz / sampleRate;
k2 = Sample(6.5451144600705975) * x +  Sample(20.46391326872472) * x * x;
"""
