import scipy.signal as signal
import soundfile
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def testAsymmetricFIR(save=False):
    fir = np.linspace(1, 0, 2**15)
    # fir /= np.sum(fir)

    if save:
        soundfile.write("data/asym_fir.wav", fir, 48000, subtype="FLOAT")

    cpp, fs = soundfile.read("snd/asym_fir_cpp.wav")

    sig = np.zeros(48000)
    sig[0] = 1
    conv = signal.convolve(sig, fir)

    plt.plot(conv, label="SciPy")
    plt.plot(cpp, label="C++")
    plt.show()

def plotAll():
    for path in Path("snd").glob("*.wav"):
        sig, fs = soundfile.read(str(path))

        plt.figure()
        plt.plot(sig, lw=1, color="red", alpha=0.75, label=path.stem)
        plt.grid()
        plt.legend()
    plt.show()

def compareFir():
    fir_real, fs = soundfile.read("snd/fir_real.wav")
    fir_imag, fs = soundfile.read("snd/fir_imag.wav")

    with open("fir.json", "r", encoding="utf-8") as fi:
        source = np.array(json.load(fi))

    spec = np.fft.rfft(source / (2 * len(source)), 2 * len(source))
    src_real = spec.real
    src_imag = spec.imag

    print(len(src_real), len(fir_real))

    plt.figure()
    plt.title("real")
    # plt.plot(fir_real, alpha=0.5, label="C++")
    # plt.plot(src_real, alpha=0.5, label="src")
    plt.plot(fir_real - src_real, alpha=0.5, label="diff")
    plt.grid()
    plt.legend()

    plt.figure()
    plt.title("imag")
    # plt.plot(fir_imag, alpha=0.5, label="C++")
    # plt.plot(src_imag, alpha=0.5, label="src")
    plt.plot(fir_imag - src_imag, alpha=0.5, label="diff")
    plt.grid()
    plt.legend()

    plt.show()

def firlp2hp(lp):
    """
    `index` search is heuristic. Check if something is wrong.
    """
    hp = -lp
    index = np.argmax(np.abs(hp))
    hp[index] += 1
    return hp

def compareCoefficient():
    coef, fs = soundfile.read("snd/coefficient.wav")

    target = signal.firwin(2047, 1000, window="nuttall", fs=48000)
    target = firlp2hp(target) / len(coef)

    plt.plot(coef, lw=1, color="blue", alpha=0.75, label="coef")
    plt.plot(target, lw=1, color="red", alpha=0.75, label="target")
    # plt.plot(target - coef[:-1], lw=1, color="red", alpha=0.75, label="diff")
    plt.grid()
    plt.legend()
    plt.show()

testAsymmetricFIR()
# plotAll()
# compareFir()
# compareCoefficient()
