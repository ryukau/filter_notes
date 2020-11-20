"""
パラメータを変えたいくつかの SOCP FIR が EBU TECH 3341 のテストを通るか調べるコード。
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile
import json
from pathlib import Path

from fractionaldelaysocp import *

def applySocpFir(sig, socpFirList):
    """FIR fractional delay filter designed by solving second order cone problem."""
    out = []
    for fir in socpFirList:
        out.append(signal.convolve(sig, fir, mode="same"))
    out = np.ravel(np.array(out).T)
    delay = len(socpFirList[0]) // 2
    index = np.linspace(-1, len(sig) - 1, len(out), endpoint=False)
    return out, index

def test(name, table):
    datapaths = [
        "data/ebu-loudness-test-setv05/seq-3341-15-24bit.wav.wav",
        "data/ebu-loudness-test-setv05/seq-3341-16-24bit.wav.wav",
        "data/ebu-loudness-test-setv05/seq-3341-17-24bit.wav.wav",
        "data/ebu-loudness-test-setv05/seq-3341-18-24bit.wav.wav",
        "data/ebu-loudness-test-setv05/seq-3341-19-24bit.wav.wav",
        "data/ebu-loudness-test-setv05/seq-3341-20-24bit.wav.wav",
        "data/ebu-loudness-test-setv05/seq-3341-21-24bit.wav.wav",
        "data/ebu-loudness-test-setv05/seq-3341-22-24bit.wav.wav",
        "data/ebu-loudness-test-setv05/seq-3341-23-24bit.wav.wav",
    ]

    targets = [
        (-6.0, +0.2, -0.4),  # 15
        (-6.0, +0.2, -0.4),  # 16
        (-6.0, +0.2, -0.4),  # 17
        (-6.0, +0.2, -0.4),  # 18
        (+3.0, +0.2, -0.4),  # 19
        (0.0, +0.2, -0.4),  # 20
        (0.0, +0.2, -0.4),  # 21
        (0.0, +0.2, -0.4),  # 22
        (0.0, +0.2, -0.4),  # 23
    ]

    isPassed = True
    for path, target in zip(datapaths, targets):
        data, samplerate = soundfile.read(path, always_2d=True)
        data = data.T[0]

        out, _ = applySocpFir(data, table)
        truepeak = np.max(np.abs(out))
        # truepeak = max(np.max(np.abs(out)), np.max(np.abs(data)))
        dbtp = 20 * np.log10(truepeak)

        targetMax = target[0] + target[1]
        targetMin = target[0] + target[2]
        status = dbtp <= targetMax and dbtp >= targetMin

        if status == False:
            print(name, Path(path).stem, dbtp, status)
            isPassed = False
            break
    return isPassed

def testRange(arrLength, arrOmegaMax):
    oversample = 4

    data = {}
    for length in arrLength:
        statusList = {}
        for omega_max in arrOmegaMax:
            name = f"FIR Length={length}, omega_max={omega_max:.3f}"
            table = createTable(length, oversample + 1, omega_max=omega_max)["table"]
            statusList[f"{omega_max:.3f}"] = test(name, table)
        data[str(length)] = statusList
    with open("socp_ebutest_status.json", "w", encoding="utf-8") as fi:
        json.dump(data, fi)
    return data

def convertToScatterData(data):
    dataX = []
    dataY = []
    for i0, vv in enumerate(data.values()):
        for i1, value in enumerate(vv.values()):
            if value:
                dataX.append(i0)
                dataY.append(i1)
    return (dataX, dataY)

def plotTestStatus(testStatus, xTick, yTick):
    dataX, dataY = convertToScatterData(testStatus)

    plt.figure(figsize=(4.5, 5), tight_layout=True)
    plt.title("Test Result, 4x Oversampling")
    plt.scatter(dataX, dataY, marker="o", s=96, zorder=2)
    plt.xticks(range(len(xTick)), xTick)
    plt.yticks(range(len(yTick)), yTick)
    plt.xlabel("FIR Length")
    plt.ylabel("omega_max")
    plt.gca().xaxis.tick_bottom()
    plt.grid()

def getPlotData(key, errors, testStatus):
    bs1770error = {
        "truepeak": 0.0011130596765247989,
        "overread": 0.0005132591093803915,
        "underread": -0.0005998005671444072
    }

    bs1770errorKey = abs(bs1770error[key])

    data = {}
    for length in arrLength:
        yList = {}
        strLength = str(length)
        for omega_max in arrOmegaMax:
            strOmegaMax = f"{omega_max:.3f}"
            yList[strOmegaMax] = (
                testStatus[strLength][strOmegaMax] and
                abs(errors[strLength][strOmegaMax][key]) < bs1770errorKey)
        data[str(length)] = yList
    return data

def plotAll(testStatus, dataTruepeak, dataOverread, dataUnderread):
    cmap = plt.get_cmap("viridis")
    plt.figure(figsize=(7, 8), tight_layout=True)
    plt.title("Finding SOCP FIR Parameter, 4x Oversampling, Marked is better.")
    plt.scatter(
        *convertToScatterData(testStatus),
        label="EBU TECH 3341 Test",
        marker=7,
        s=192,
        zorder=2,
        color=cmap(0.0),
    )
    plt.scatter(
        *convertToScatterData(dataTruepeak),
        label="Truepeak",
        marker=6,
        s=192,
        zorder=2,
        color=cmap(0.3),
    )
    plt.scatter(
        *convertToScatterData(dataOverread),
        label="Overread",
        marker=4,
        s=192,
        zorder=2,
        color=cmap(0.6),
    )
    plt.scatter(
        *convertToScatterData(dataUnderread),
        label="Underread",
        marker=5,
        s=192,
        zorder=2,
        color=cmap(0.9),
    )
    plt.xticks(range(len(arrLength)), arrLength)
    plt.yticks(range(len(tickOmegaMax)), tickOmegaMax)
    plt.xlabel("FIR Length")
    plt.ylabel("omega_max")
    plt.grid()
    plt.legend()
    plt.show()

def plotReasonable(testStatus, dataTruepeak, dataOverread, dataUnderread):
    dataX = []
    dataY = []
    for i0, (length, vv) in enumerate(testStatus.items()):
        for i1, (omega_max, value) in enumerate(vv.items()):
            if (value and dataTruepeak[length][omega_max] and
                    dataOverread[length][omega_max] and dataUnderread[length][omega_max]):
                dataX.append(i0)
                dataY.append(i1)
    plt.figure(figsize=(4.2, 5), tight_layout=True)
    plt.title("Reasonable SOCP FIR Parameter,\n4x Oversampling")
    plt.scatter(dataX, dataY, marker="o", s=64, zorder=2)
    plt.xticks(range(len(arrLength)), arrLength)
    plt.yticks(range(len(tickOmegaMax)), tickOmegaMax)
    plt.xlabel("FIR Length")
    plt.ylabel("omega_max")
    plt.grid()
    plt.show()

if __name__ == "__main__":
    arrLength = np.arange(3, 17)
    arrOmegaMax = np.arange(0.5, 0.91, 0.025)
    tickOmegaMax = [f"{x:.3f}" for x in arrOmegaMax]

    testStatus = testRange(arrLength, arrOmegaMax)
    # plotTestStatus(testStatus, arrLength, tickOmegaMax)

    with open("socp_ebutest_status.json", "r", encoding="utf-8") as fi:
        testStatus = json.load(fi)
    with open("error_socp_omega_max.json", "r", encoding="utf-8") as fi:
        errors = json.load(fi)

    dataTruepeak = getPlotData("truepeak", errors, testStatus)
    dataOverread = getPlotData("overread", errors, testStatus)
    dataUnderread = getPlotData("underread", errors, testStatus)

    plotAll(testStatus, dataTruepeak, dataOverread, dataUnderread)
    plotReasonable(testStatus, dataTruepeak, dataOverread, dataUnderread)
