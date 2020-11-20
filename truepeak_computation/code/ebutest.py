import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import scipy.special as special
import soundfile
import json
from pathlib import Path

class LagrangeFD:
    def __init__(self, order):
        self.order = order
        self.xd = np.zeros(order)
        self.diff = np.zeros(order)

    def reset(self):
        self.xd = 0
        self.diff = 0

    def push(self, x0):
        self.diff[0] = x0 - self.xd[0]
        self.xd[0] = x0
        for i in range(1, self.order):
            self.diff[i] = self.diff[i - 1] - self.xd[i]
            self.xd[i] = self.diff[i - 1]

    def at(self, fraction):
        delta = fraction + (self.order - 1) / 2
        sig = 0

        for i in range(self.order, 0, -1):
            sig = ((i - 1) - delta) / i * (self.diff[i - 1] + sig)

        return sig + self.xd[0]

    def process(self, x0, fraction):
        self.push(x0)
        return at(fraction)

def thiranAllpass(order, fraction):
    """
    Return transfer function coefficients (b, a).
    fraction > 0.
    Delay is order + fraction.
    """
    a = np.empty(order)
    for k in range(order):
        a[k] = (-1)**k * special.binom(order, k)
        for n in range(order):
            a[k] *= (fraction + n) / (fraction + k + n)
    return (a[::-1], a)

def iturBs1770Fir():
    return [
        [
            0.001708984375, 0.010986328125, -0.0196533203125, 0.033203125,
            -0.0594482421875, 0.1373291015625, 0.97216796875, -0.102294921875,
            0.047607421875, -0.026611328125, 0.014892578125, -0.00830078125
        ],
        [
            -0.0291748046875, 0.029296875, -0.0517578125, 0.089111328125, -0.16650390625,
            0.465087890625, 0.77978515625, -0.2003173828125, 0.1015625, -0.0582275390625,
            0.0330810546875, -0.0189208984375
        ],
        [
            -0.0189208984375, 0.0330810546875, -0.0582275390625, 0.1015625,
            -0.2003173828125, 0.77978515625, 0.465087890625, -0.16650390625,
            0.089111328125, -0.0517578125, 0.029296875, -0.0291748046875
        ],
        [
            -0.00830078125, 0.014892578125, -0.026611328125, 0.047607421875,
            -0.102294921875, 0.97216796875, 0.1373291015625, -0.0594482421875,
            0.033203125, -0.0196533203125, 0.010986328125, 0.001708984375
        ],
    ]

def socpFir():
    with open("socp.json", "r", encoding="utf-8") as fi:
        data = json.load(fi)
    return list(reversed(data["table"][1:]))

def applyFullSinc(sig, oversample=32):
    out = []
    index = np.arange(-len(sig), len(sig))
    for fraction in np.linspace(0, 1, oversample, endpoint=False):
        sinc = np.sinc(index + fraction)
        out.append(signal.convolve(sig, sinc, mode="same"))
    out = np.ravel(np.array(out).T)
    index = np.linspace(-1, len(sig) - 1, len(out), endpoint=False)
    return out, index

def applyPartialSinc(sig, half=9, oversample=32):
    out = []
    index = np.arange(-half, half + 1)
    for fraction in np.linspace(0, 1, oversample, endpoint=False):
        sinc = np.sinc(index + fraction)
        out.append(signal.convolve(sig, sinc, mode="same"))
    out = np.ravel(np.array(out).T)
    index = np.linspace(-1, len(sig) - 1, len(out), endpoint=False)
    return out, index

def applyLagrangeInterpolation(sig, order=11):
    fd = LagrangeFD(order)
    order1 = order + 1
    out = np.empty(len(sig) * order1)
    fraction = np.linspace(1, 0, order1, endpoint=False)
    idx = 0
    for x in sig:
        fd.push(x)
        for frac in fraction:
            out[idx] = fd.at(frac)
            idx += 1
    delay = order1 // 2
    index = np.linspace(-delay, len(sig) - delay, len(out), endpoint=False)
    return out, index

def applyBs1770Fir(sig):
    """FIR from ITU-R BS.1770-4, p.17 (Annex 2-3)."""
    out = []
    for fir in iturBs1770Fir():
        out.append(signal.convolve(sig, fir, mode="same"))
    out = np.ravel(np.array(out).T)
    index = np.linspace(-1, len(sig) - 1, len(out), endpoint=False)
    return out, index

def applySocpFir(sig):
    """FIR fractional delay filter designed by solving second order cone problem."""
    out = []
    socpFirList = socpFir()
    for fir in socpFirList:
        out.append(signal.convolve(sig, fir, mode="same"))
    out = np.ravel(np.array(out).T)
    delay = len(socpFirList[0]) // 2
    index = np.linspace(-1, len(sig) - 1, len(out), endpoint=False)
    return out, index

def applyBessel(sig, order=8):
    delays = np.linspace(0, 1, 16, endpoint=False) + 1
    out = []
    for dly in delays:
        sos = signal.bessel(order, 2 / (np.pi * dly), norm="delay", output="sos")
        out.append(signal.sosfilt(sos, sig))
    out = np.ravel(np.array(out).T)
    index = np.linspace(0, len(sig) - 1, len(out), endpoint=False)
    return out, index

def applyThiran(sig, order=12, oversample=8):
    out = []
    for fraction in np.linspace(0, 1, oversample, endpoint=False)[1:]:
        b, a = thiranAllpass(order, fraction)
        sos = signal.tf2sos(b, a)
        out.append(signal.sosfilt(sos, sig))
    delay = order // 2
    index = np.linspace(-delay, len(sig) - delay, len(out), endpoint=False)
    return out, index

def test(name, filterFunc):
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
        # "data/worstsinc/worst_48000Hz_01sec.wav",
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
        # (17.696852601831182, +0.2, -0.4),  # worstsinc 1 sec
    ]

    print(f"--- {name}")
    for path, target in zip(datapaths, targets):
        data, samplerate = soundfile.read(path, always_2d=True)
        data = data.T[0]

        out, _ = filterFunc(data)
        dataAbs = np.abs(data)
        truepeak = np.max(np.abs(out))
        dbtp = 20 * np.log10(truepeak)

        targetMax = target[0] + target[1]
        targetMin = target[0] + target[2]
        status = dbtp <= targetMax and dbtp >= targetMin

        print(Path(path).stem, dbtp, status)
    print()

for func in [
    ("Full Sinc", applyFullSinc),
    ("Partial Sinc", applyPartialSinc),
    ("Lagrange", applyLagrangeInterpolation),
    ("BS.1770 FIR", applyBs1770Fir),
    ("SOCP FIR", applySocpFir),
    ("Bessel", applyBessel),
    ("Thiran AP", applyThiran),
]:
    test(*func)
