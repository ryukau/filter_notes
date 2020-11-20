import matplotlib.pyplot as pyplot
import numpy as np
import scipy.signal as signal
import scipy.special as special
import soundfile

def thiran(order, fraction):
    """
    Return transfer function coefficients (a_reversed, a).
    fraction in [0, 1].
    """
    a = np.empty(order)
    for k in range(order):
        a[k] = (-1)**k * special.binom(order, k)
        for n in range(order):
            a[k] *= (fraction - order + n) / (fraction - order + k + n)

    return (a[::-1], a)

def plot_thiran_phase(order=4):
    fraction = np.linspace(order, order + 1, 16)
    print(fraction[5])
    cmap = pyplot.get_cmap("viridis")
    for idx, frac in enumerate(fraction):
        coef = thiran(order, frac)
        freq, gd = signal.group_delay(coef, w=1024)
        color = idx / len(fraction)
        if idx == 0 or idx == 15:
            pyplot.plot(freq, gd, color=cmap(color), alpha=0.5, label=f"{idx}")
        else:
            pyplot.plot(freq, gd, color=cmap(color), alpha=0.5)
    pyplot.xlabel("Frequency [rad/sample]")
    pyplot.ylabel("Group Delay [samples]")
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

def testApply():
    samplerate = 44100
    duration = 1  # seconds
    frequency = 100

    time = np.linspace(0, duration, duration * samplerate)
    phase = (frequency * time)

    saw = signal.sawtooth(phase, 1.0)

    coef = thiran(16, 16.5)
    print(coef)
    sos = signal.tf2sos(*coef)
    sig = signal.sosfilt(sos, saw)
    print(f"Converge?: {np.all(np.isfinite(sig))}")
    soundfile.write("out.wav", sig, samplerate, subtype="FLOAT")

def printSos(order=12, oversample=5):
    for fraction in np.linspace(0, 1, oversample, endpoint=False)[1:]:
        print(fraction)
        b, a = thiran(order, order + fraction)
        sos = signal.tf2sos(b, a)
        print("{")
        for section in sos:
            text = "  {"

            for idx, value in enumerate(section):
                if idx == 3:
                    continue
                if idx != len(section) - 1:
                    text += f"Sample({value}), "
                else:
                    text += f"Sample({value})"
            print(text + "},")
        print("}")

printSos()
# testApply()
# plot_thiran_phase()
