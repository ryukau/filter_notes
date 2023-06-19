import json
import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def toDecibel(response):
    return 20 * np.log10(np.maximum(np.abs(response), 1e-7))


def plotResponse(key, bands):
    ir_sum = np.empty_like(bands[0])
    for idx, bd in enumerate(bands):
        ir_sum += np.array(bd)
        freq, resp = signal.freqz(bd, fs=2)
        plt.plot(freq, toDecibel(resp), alpha=0.5, label=f"band{idx}")

    freq, resp = signal.freqz(ir_sum, fs=2)
    plt.plot(freq, toDecibel(resp), alpha=0.5, label="sum")

    plt.axvline(1000 / 24000, color="black", ls="--")
    plt.title(key)
    plt.xscale("log")
    plt.ylim([-60, 10])
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


def plotImpulse(key, bands):
    ir_sum = np.zeros_like(bands[0])
    for idx, bd in enumerate(bands):
        ir_sum += np.array(bd)
        plt.plot(bd, alpha=0.5, label=f"band{idx}")
    # plt.plot(ir_sum, alpha=0.5, label="sum")
    plt.title(key)
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    with open("ir.json", "r", encoding="utf-8") as fi:
        data = json.load(fi)

    for key, bands in data.items():
        plotResponse(key, bands)
        # plotImpulse(key, bands)
