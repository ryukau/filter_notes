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
    args = (path, tables)
    tables = {
        "firName": [
            [b00, b01, ...], # FIR 0
            [b10, b11, ...], # FIR 1
            ...
        ],
        ...
    }
    """
    path = args[0]

    if not path.is_file():
        return None

    try:
        data, samplerate = soundfile.read(str(path), always_2d=True)
    except:
        print(f"Skipping: {path}")
        return None
    print(f"Processing: {path}")

    result = {}
    data = data.T[0]
    for name, table in args[1].items():
        truepeak = applySocpFir(data, table)
        dbtp = 20 * np.log10(truepeak)
        result[name] = {"True-peak": truepeak, "dBTP": dbtp}

    return (str(path.as_posix()), result)

def dumpTable():
    tables = {}
    omega_max = 0.5
    for idx in range(4, 17):
        tbl = createTable(12, idx + 1, omega_max=omega_max)
        tables[f"SocpOversampleEven{idx:02d}"] = tbl["table"].tolist()[:-1]
    for idx in range(4, 17):
        tbl = createTable(13, idx + 1, omega_max=omega_max)
        tables[f"SocpOversampleOdd{idx:02d}"] = tbl["table"].tolist()[:-1]
    for idx in range(12, 33):
        tbl = createTable(idx, 6, omega_max=omega_max)
        tables[f"SocpLengthOdd{idx:02d}"] = tbl["table"].tolist()[:-1]
    for idx in range(12, 33):
        tbl = createTable(idx, 7, omega_max=omega_max)
        tables[f"SocpLengthEven{idx:02d}"] = tbl["table"].tolist()[:-1]

    with open("measure_socp_filter.json", "w", encoding="utf-8") as fi:
        json.dump(tables, fi)

if __name__ == "__main__":
    dumpTable()

    with open("measure_socp_filter.json", "r", encoding="utf-8") as fi:
        tables = json.load(fi)

    result = {}
    paths = [(p, tables) for p in Path("data").glob("**/*")]
    with Pool() as pool:
        for data in pool.imap_unordered(job, paths):
            if data is None:
                continue
            result[data[0]] = data[1]

    with open("measure_socp.json", "w", encoding="utf-8") as fi:
        json.dump(result, fi)
