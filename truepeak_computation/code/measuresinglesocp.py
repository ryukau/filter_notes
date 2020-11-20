"""
SOCP FIR の誤差を計算するコード。
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile
import json
from multiprocessing import Pool
from pathlib import Path

from fractionaldelaysocp import *

def applySocpFir(sig, table):
    """FIR fractional delay filter designed by solving second order cone problem."""
    out = [np.max(np.abs(sig))]
    for fir in table:
        socp = signal.convolve(sig, fir, mode="same")
        out.append(np.max(np.abs(socp)))
    return np.max(out)

def job(args):
    """
    args = (path, table)
    tables = [
        [b00, b01, ...], # FIR 0
        [b10, b11, ...], # FIR 1
        ...
    ]
    """
    path = args[0]

    if not path.is_file():
        return None

    try:
        data, samplerate = soundfile.read(str(path), always_2d=True)
    except:
        # print(f"Skipping: {path}")
        return None
    # print(f"Processing: {path}")

    data = data.T[0]
    truepeak = applySocpFir(data, args[1])
    dbtp = 20 * np.log10(truepeak)
    result = {"True-peak": truepeak, "dBTP": dbtp}

    return (str(path.as_posix()), result)

def getTruepeak(data):
    return (data["True-peak"], data["dBTP"])

def measureSocp(length=7, oversample=4, omega_max=0.65):
    table = createTable(length, oversample + 1, omega_max=omega_max)["table"][:-1]

    dataSocp = {}
    paths = [(p, table) for p in Path("data").glob("**/*")]
    with Pool() as pool:
        for data in pool.imap_unordered(job, paths):
            if data is None:
                continue
            dataSocp[data[0]] = data[1]

    with open("measure_single_socp.json", "w", encoding="utf-8") as fi:
        json.dump(dataSocp, fi)

    with open("measure_sinc.json") as fi:
        dataSinc = json.load(fi)

    truepeak = 0
    overread = 0
    underread = 0
    for dataName in dataSinc.keys():
        if dataName == "data/worstsinc/worst_48000Hz_01sec.wav":
            continue
        truepeakSinc, _ = getTruepeak(dataSinc[dataName])
        truepeakSocp, _ = getTruepeak(dataSocp[dataName])

        diffTp = truepeakSocp - truepeakSinc
        truepeak += abs(diffTp)
        if diffTp > 0:
            overread += diffTp
        elif diffTp < 0:
            underread += diffTp
    return {
        "truepeak": truepeak / len(dataSinc),
        "overread": overread / len(dataSinc),
        "underread": underread / len(dataSinc),
    }

def measureSocpVaryingOmegaMax():
    oversample = 4

    result = {}
    for length in np.arange(3, 17):
        resultLength = {}
        for omega_max in np.arange(0.5, 0.91, 0.025):
            print(f"Processing: {length}, {omega_max}")
            resultLength[f"{omega_max:.3f}"] = measureSocp(length, oversample, omega_max)
        result[str(length)] = resultLength

    with open("error_socp_omega_max.json", "w", encoding="utf-8") as fi:
        json.dump(result, fi)

def plot(keyName: str):
    with open("error_socp_omega_max.json", "r", encoding="utf-8") as fi:
        data = json.load(fi)

    nPlot = len(data.keys())
    nRows = nPlot // 2
    nCols = nPlot // nRows

    plt.figure(1, figsize=(9, 3 * nRows), tight_layout=True)
    plt.suptitle(f"{keyName.capitalize()} Error of SOCP FIR")
    for idx, (length, dataOmegaMax) in enumerate(data.items()):
        plt.subplot(nRows, nCols, idx + 1)
        xx = []
        yy = []
        for omega_max, values in dataOmegaMax.items():
            if float(omega_max) > 0.81:
                continue
            xx.append(f"{float(omega_max):.2f}")
            yy.append(values[keyName])
        plt.title(f"FIR Length {length}")
        plt.bar(xx, yy, color="red")
        if idx % nCols == 0:
            plt.ylabel("Error")
        if idx // nCols == nRows - 1:
            plt.xlabel("omega_max")
    plt.show()

if __name__ == "__main__":
    measureSocpVaryingOmegaMax()
    # plot("truepeak")
    # plot("overread")
    # plot("underread")
