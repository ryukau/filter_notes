"""
有限の長さの信号について、最悪の場合のトゥルーピークの大きさを計算するコード。
"""

import numpy as np
from scipy.special import psi

def calcWorstTruePeak(nSample, fraction=0.5):
    m = nSample // 2
    t = fraction
    A = (psi(t + m + 1) - psi(t)) * np.sin(np.pi * t) / np.pi
    B = (psi(-t + m + 2) - psi(1 - t)) * np.sin(np.pi * (t - 1)) / np.pi
    return A - B

fs = 48000

print(calcWorstTruePeak(fs))
print(calcWorstTruePeak(60 * fs))
print(calcWorstTruePeak(60 * 60 * fs))
print(calcWorstTruePeak(24 * 60 * 60 * fs))
