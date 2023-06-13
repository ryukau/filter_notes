import json
import matplotlib.pyplot as plt
import numpy as np


def lpCut(freq):
    y = 1 - np.cos(2 * np.pi * freq)
    return np.sqrt((y + 2) * y) - y


def apCut(freq):
    t = np.tan(np.pi * freq)
    return (t - 1) / (t + 1)


def plot(data):
    freq = np.array(data["normalizedFreq"])
    reso = np.array(data["resonance"])

    c1 = lpCut(freq)
    c2 = apCut(freq)

    # prod = 0.5 * np.max(reso) * ((1 - c1) * c2 + c1 + 1)

    plt.plot(freq, reso, label="Actual")
    # plt.plot(freq, prod, label="prod")

    # plt.plot(np.diff(reso, 2), label="Actual")
    # plt.plot(np.diff(prod, 2), label="prod")

    plt.legend()
    plt.grid()
    plt.show()


def reduce(data):
    freq = np.array(data["normalizedFreq"])
    reso = np.array(data["resonance"])

    nReduced = 255
    step = int(np.ceil(len(freq) / nReduced))
    reducedReso = []
    reducedFreq = []
    for index in range(len(freq) - 1, 0, -step):
        reducedFreq.insert(0, freq[index])
        reducedReso.insert(0, reso[index])
    reducedFreq.insert(0, 0)
    reducedReso.insert(0, 0)

    linterp = np.interp(freq, reducedFreq, reducedReso)
    fix = np.max(linterp[1:] / reso[1:]) + np.finfo(np.float64).eps
    linterp /= fix

    if np.any(linterp > reso):
        print("Unstable")
        exit()

    # plt.plot(freq, reso, label="Actual")
    # plt.plot(freq, linterp, label="Reduced")
    # plt.legend()
    # plt.grid()
    # plt.show()

    return (reducedFreq, reducedReso)


def writeTable(resonance):
    text = ""

    text += f"""// This file is generated.
#pragma once
#include <array>

namespace ResonantEmaUtil {{
template<typename T> struct ResonanceTable {{
    static constexpr std::array<T, {len(resonance)}> qTable{{"""

    for value in resonance:
        text += f"T({value}),"

    text += f"""}};
}};
}} // namespace ResonantEmaUtil
"""
    with open("qtable.hpp", "w", encoding="utf-8") as fi:
        fi.write(text)


def processResult():
    with open("table_linear.json", "r", encoding="utf-8") as fi:
        data = json.load(fi)

    # plot(data)
    writeTable(data["resonance"])


def compareResults():
    with open("table.json", "r", encoding="utf-8") as fi:
        dataLinRegress = json.load(fi)
    with open("result2.json", "r", encoding="utf-8") as fi:
        dataPeak = json.load(fi)

    freqL = dataLinRegress["normalizedFreq"]
    freqP = dataPeak["normalizedFreq"]
    if np.any(freqL != freqP):
        print("Discrepancy in normalized frequencies.")
        exit()

    diff = np.abs(resoL - resoP)
    cutoff = int(len(freqL) * 100 / 48000)

    freq = freqL[:cutoff]
    resoL = np.array(dataLinRegress["resonance"][:cutoff])
    resoP = np.array(dataPeak["resonance"][:cutoff])

    # plt.plot(freq, resoL, label="linregress")
    # plt.plot(freq, resoP, label="peak")
    plt.plot(freq, resoL - resoP, label="diff")
    # plt.yscale("log")
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    processResult()
    # compareResults()
