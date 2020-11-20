"""
Plot responses of FIR filter from ITU-R BS.1770-4, p17 (Annex 2-3).
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import json

def plot(firTable, ylim):
    fig, ax = plt.subplots(3, 1)
    fig.set_size_inches(6, 8)
    fig.set_tight_layout(True)

    ax[0].set_title("Amplitude Response")
    ax[0].set_ylabel("Amplitude [dB]")
    ax[1].set_title("Phase Response")
    ax[1].set_ylabel("Phase [degree]")
    ax[2].set_title("Group Delay")
    ax[2].set_ylabel("Group Delay [sample]")
    ax[2].set_xlabel("Frequency [rad/sample]")
    ax[2].set_ylim(ylim)

    cmap = plt.get_cmap("plasma")
    for index, fir in enumerate(firTable):
        freq, mag, phase = signal.dbode((fir, 1, 1), n=1024)
        freq, gd = signal.group_delay((fir, 1), w=1024)

        cmapValue = index / len(firTable)
        ax[0].plot(freq, mag, alpha=0.75, color=cmap(cmapValue), label=f"{index}")
        ax[1].plot(freq, phase, alpha=0.75, color=cmap(cmapValue), label=f"{index}")
        ax[2].plot(freq, gd, alpha=0.75, color=cmap(cmapValue), label=f"{index}")

    for axis in ax:
        axis.grid()
        axis.legend()
    plt.show()

def plotBs1770fir():
    with open("bs1770fir.json", "r", encoding="utf-8") as fi:
        data = json.load(fi)

    bs1770fir = []
    for key, value in data.items():
        bs1770fir.append(value)
    plot(bs1770fir, ylim=(3.5, 7.5))

def plotSocpfir():
    with open("socp.json", "r", encoding="utf-8") as fi:
        data = json.load(fi)
    yMid = len(data["table"][0]) / 2 - 0.5
    plot(data["table"], ylim=(yMid - 2, yMid + 2))

# plotBs1770fir()
plotSocpfir()
