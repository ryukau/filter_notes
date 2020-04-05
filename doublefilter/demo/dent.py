import numpy as np
import scipy.signal as signal
import scipy.optimize as optimize
import matplotlib.pyplot as pyplot

def shelveK1(k1, k2):
    k1Gain = 0.69
    k1Delta = 0.31  # 1 - k1Gain.
    B_2 = 0.63  # Tuning boundary (xi_{k_2}).
    B_3 = 0.635

    if k2 < B_2:
        k1 *= k1Gain
    elif k2 >= B_2 and k2 < B_3:
        k1 *= k1Gain + k1Delta * (k2 - B_2) / (B_3 - B_2)
    return k1

def dentK1(k1, k2):
    k1Gain = 0.69
    k1Delta = 0.31  # 1 - k1Gain.
    B_0 = 0.61
    B_1 = 0.625
    B_2 = 0.63  # Tuning boundary (xi_{k_2}).
    B_3 = 0.635

    if k2 >= B_0 and k2 < B_1:
        k1 *= 1 - k1Delta * (k2 - B_0) / (B_1 - B_0)
    elif k2 >= B_1 and k2 < B_2:
        k1 *= k1Gain
    elif k2 >= B_2 and k2 < B_3:
        k1 *= k1Gain + k1Delta * (k2 - B_2) / (B_3 - B_2)

    return k1 * 0.7

def resonanceCurve(resonance, k2, tuningFunc):
    C0 = -0.0049691265927442885
    C1 = -471.738128187657
    C2 = 1432.5662635997667
    C3 = 345.2853784111966
    C4 = -4454.40786711102
    C5 = 3468.062963176107
    k1 = np.where(
        k2 < 0.6295160864148501, resonance * np.pi,
        resonance *
        (C0 + 1 /
         (C1 + C2 * k2 + C3 * k2 * k2 + C4 * k2 * k2 * k2 + C5 * k2 * k2 * k2 * k2)))
    return (k1, np.array([tuningFunc(k1v, k2v) for k1v, k2v in zip(k1, k2)]))

def plotK1Curve(tuningFunc):
    reso = 1
    k2 = np.linspace(0, 1, 4096)

    k1Raw, k1Shelved = resonanceCurve(reso, k2, tuningFunc)

    pyplot.figure(figsize=(6, 3))
    pyplot.title(f"Resonance={reso}")
    pyplot.plot(k1Raw, label="Raw")
    pyplot.plot(k1Shelved, label="Shelved")
    pyplot.xlabel("k2")
    pyplot.ylabel("k1")
    pyplot.grid()
    pyplot.legend()
    pyplot.tight_layout()
    pyplot.show()

if __name__ == "__main__":
    plotK1Curve(shelveK1)
    plotK1Curve(dentK1)
