import numpy as np
import scipy.signal as signal
import scipy.optimize as optimize
import matplotlib.pyplot as pyplot
import json

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

def plotPole(transferFunction, k1, k2):
    cmap = pyplot.get_cmap("viridis")
    b, a = transferFunction(k1, k2)
    z, p, k = signal.tf2zpk(b, a)
    z *= k
    p *= k
    if np.max(np.abs(p)) >= 1:
        print("Unstable")
    circle1 = pyplot.Circle((0, 0), 1, fill=False, color="black", zorder=2)
    pyplot.gcf().gca().add_artist(circle1)
    pyplot.plot(np.real(z), np.imag(z), "xb", zorder=3, label="Zeros")
    pyplot.plot(np.real(p), np.imag(p), "or", zorder=3, label="Poles")
    pyplot.legend(loc=2)
    pyplot.title(f"Pole / Zero Plot k={k:.3f}")
    pyplot.ylabel("Real")
    pyplot.xlabel("Imaginary")
    pyplot.xlim((-2, 2))
    pyplot.ylim((-2, 2))
    pyplot.grid(zorder=1)
    pyplot.tight_layout()
    pyplot.show()

def stabilityPlot(transferFunction):
    pyplot.figure(figsize=(6, 3))
    maxIteration = 1024
    nK2 = 256
    data = []
    for idx, k2 in enumerate(np.linspace(0.1, 1, nK2)):
        k1 = 65536
        delta = k1 / 2
        jdx = 0
        while jdx < maxIteration:
            b, a = transferFunction(k1, k2)
            _, pole, gain = signal.tf2zpk(b, a)
            if np.max(np.abs(pole * gain)) >= 1:
                k1 -= delta
            else:
                k1 += delta
            jdx += 1
            delta *= 0.5
        data.append((k1, k2))
        print(idx, k1, k2)
    k1, k2 = zip(*data)
    pyplot.scatter(k2, k1, s=4, zorder=2)
    k1Array, k2Array = zip(*data)
    # with open("k1_k2.json", "w") as fi:
    #     json.dump({"k1": k1Array, "k2": k2Array}, fi, indent=2)
    pyplot.title("k1-k2 Plot")
    pyplot.ylabel("k1")
    pyplot.xlabel("k2")
    pyplot.grid(zorder=1)
    pyplot.tight_layout()
    pyplot.show()

def fitK1K2():
    pyplot.figure(figsize=(6, 3))
    with open("k1_k2_0.63.json", "r") as fi:
        data = json.load(fi)
    k1 = np.array(data["k1"])
    k2 = np.array(data["k2"])

    # pyplot.plot(k2, 1 / k2, lw=1, alpha=0.5, color="red", label="fit")

    def curve_func(x, a0, a1, a2, a3, a4):
        return 1 / (a0 + a1 * x + a2 * x * x + a3 * x * x * x + a4 * x * x * x * x)

    popt, _ = optimize.curve_fit(curve_func, np.copy(k2), np.copy(k1))
    print([p for p in popt])
    plotX = np.linspace(k2[0], k2[-1], 2**10)
    curve = curve_func(plotX, *popt) - 0.01965707063273614
    pyplot.plot(plotX, curve, lw=1, ls="--", alpha=0.75, color="red", label="fit")

    # error = (curve_func(k2, *popt) - 0.01965707063273614 - k1)
    # print(f"Max Positive Error: {np.max(error)} at {np.argmax(error)}")
    # print(curve[0], curve[-1])

    # print(curve_func(0.2 * np.pi, *popt))

    pyplot.plot(k2, k1, lw=1, alpha=0.75, color="blue", label="Real")
    pyplot.ylabel("k1")
    pyplot.xlabel("k2")
    # pyplot.ylim((-10, 110))
    pyplot.grid(zorder=1)
    pyplot.legend()
    pyplot.tight_layout()
    pyplot.show()

# plotPole(transferFunctionU2, 10000, 1e-05)
stabilityPlot(transferFunctionU2)
# fitK1K2()
"""
```maxima
a0: -471.738128187657;
a1: 1432.5662635997667;
a2: 345.2853784111966;
a3: -4454.40786711102;
a4: 3468.062963176107;
f(x) := 1 / (a0 + a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4);
numer: true;
solve(%pi = f(x) - 0.003127261304194378, x);
```

When x is 0.6295160864148501, f(x) - 0.0049691265927442885 is π.
"""
