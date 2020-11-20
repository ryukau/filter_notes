"""
Sinc 補間によるトゥルーピークを計算するコード。

32 bit 環境ではメモリ不足になるので注意。
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile
import json
from multiprocessing import Pool
from pathlib import Path

def getSincTruePeak(sig, oversample=32):
    out = []
    index = np.arange(-len(sig), len(sig))
    # index = np.arange(-32768, 32768)
    for fraction in np.linspace(0, 1, oversample, endpoint=False):
        sinc = np.sinc(index + fraction)
        sig = signal.convolve(sig, sinc, mode="same")
        out.append(np.max(np.abs(sig)))
    return np.max(out)

def job(path):
    if not path.is_file():
        return None

    try:
        data, samplerate = soundfile.read(str(path), always_2d=True)
    except:
        print(f"Skipping: {path}")
        return None

    data = data.T[0]
    truepeak = getSincTruePeak(data)
    dbtp = 20 * np.log10(truepeak)
    print(f"{path}: {truepeak}, {dbtp}[dBTP]")

    return (str(path.as_posix()), {"True-peak": truepeak, "dBTP": dbtp})

if __name__ == "__main__":
    result = {}
    paths = [p for p in Path("data").glob("**/*")]
    with Pool() as pool:
        for data in pool.imap_unordered(job, paths):
            if data is None:
                continue
            result[data[0]] = data[1]

    with open("measure_sinc.json", "w", encoding="utf-8") as fi:
        json.dump(result, fi)
