import json
import math
import pprint
import matplotlib.pyplot as plt

def getTruepeak(data):
    return (data["True-peak"], data["dBTP"])

def addAbsoluteError(data, result, truepeakSinc, dbtpSinc):
    truepeak, dbtp = getTruepeak(data)

    diffTp = truepeak - truepeakSinc
    result["truepeak"] += abs(diffTp)
    if diffTp > 0:
        result["overread"] += diffTp
    elif diffTp < 0:
        result["underread"] += diffTp

    result["dbtp"] += abs(dbtpSinc - dbtp)
    if "elapsedMilliSeconds" in data:
        result["elapsedMilliSeconds"] += data["ElapsedMilliSeconds"]

def emptyErrorSum():
    """
    underread: True-peak is lower than sinc interpolated value.
    overread: True-peak is higher than sinc interpolated value.
    """
    return {
        "truepeak": 0,
        "underread": 0,
        "overread": 0,
        "dbtp": 0,
        "elapsedMilliSeconds": 0,
    }

with open("measure_sinc.json") as fi:
    dataSinc = json.load(fi)
with open("measure_cpp.json") as fi:
    dataCpp = json.load(fi)
with open("measure_socp.json") as fi:
    dataSocp = json.load(fi)

dataMerged = {}
for dataName, dt in dataCpp.items():
    for method, value in dataSocp[dataName].items():
        dt[method] = value
    dataMerged[dataName] = dt

errorSum = {}

firstDataKey = "data/ebu-loudness-test-setv05/1kHz Sine -20 LUFS-16bit.wav"
for key, value in dataMerged[firstDataKey].items():
    if not isinstance(value, dict):
        continue
    errorSum[key] = emptyErrorSum()

for dataName, measured in dataMerged.items():
    if dataName == "data/worstsinc/worst_48000Hz_01sec.wav":
        continue
    truepeakSinc, dbtpSinc = getTruepeak(dataSinc[dataName])
    for method, result in errorSum.items():
        addAbsoluteError(measured[method], result, truepeakSinc, dbtpSinc)

for method, result in errorSum.items():
    result["truepeak"] /= len(dataSinc)
    result["underread"] /= len(dataSinc)
    result["overread"] /= len(dataSinc)
    result["dbtp"] /= len(dataSinc)

pp = pprint.PrettyPrinter(indent=2)
pp.pprint(errorSum)

def plotSocpResult(method, valueName, subplot):
    # plt.figure()
    plt.subplot(*subplot)
    plotName = []
    plotValue = []
    for key, value in errorSum.items():
        if method not in key:
            continue
        plotName.append(key.replace(method, ""))
        plotValue.append(value[valueName])

    minY = min(plotValue)
    maxY = max(plotValue)
    marginY = 0.05 * (maxY - minY)
    plt.ylim((minY - marginY, maxY + marginY))

    methodReplace = {
        "SocpOversampleEven":
            ("Even Length (12), Varying Oversample", "Oversampling Multiplier"),
        "SocpOversampleOdd":
            ("Odd Length (13), Varying Oversample", "Oversampling Multiplier"),
        "SocpLengthOdd": ("Odd Oversample (5), Varying Length", "FIR Length"),
        "SocpLengthEven": ("Even Oversample (6), Varying Length", "FIR Length"),
    }
    method, xlabel = methodReplace[method]

    plt.title(f"{method}, {valueName}")
    plt.bar(plotName, plotValue, color="red")
    plt.xlabel(xlabel)

    if valueName == "truepeak":
        plt.ylabel("Mean Absolute Error")
    elif valueName == "dbtp":
        plt.ylabel("Mean Absolute Error [dB]")
    else:
        plt.ylabel("Mean Error")

socpLabel = ["SocpOversampleEven", "SocpOversampleOdd", "SocpLengthEven", "SocpLengthOdd"]
valueLabel = ["truepeak", "overread", "underread", "dbtp"]

for valuelbl in valueLabel:
    plt.figure(figsize=(6, 9), tight_layout=True)
    for idx, socplbl in enumerate(socpLabel):
        plotSocpResult(socplbl, valuelbl, (len(socpLabel), 1, idx + 1))
plt.show()
