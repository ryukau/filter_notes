"""
トゥルーピークが最大になる最悪の場合の信号をレンダリングするコード。
"""

import numpy as np
import soundfile
from pathlib import Path

def check(sig):
    for i in range(1, len(sig)):
        if sig[i] == sig[i - 1]:
            print(f"Consecutive same value starts from index {i}.")
            break
    return

samplerate = 48000
fraction = 0.5  # Output signal is same if `fraction % 1 != 0`.

outDir = Path("data/worstsinc")
if not outDir.exists():
    outDir.mkdir(parents=True, exist_ok=True)

for length in (np.arange(1, 11) * samplerate):
    index = np.arange(1 - length // 2, 1 + length // 2)
    sig = np.sinc(index + fraction)
    absed = np.abs(sig)
    sig[absed < 1e-7] = 0
    sig = np.sign(sig)

    check(sig)

    path = outDir / Path(f"worst_{samplerate}Hz_{int(length/samplerate):02d}sec.wav")
    soundfile.write(str(path), sig, samplerate, subtype="FLOAT")
