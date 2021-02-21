import gc
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile
from pathlib import Path

def plot(sig, samplerate, name):
    freq, time, spectre = signal.spectrogram(
        sig,
        samplerate,
        # ("tukey", 0.33),
        "blackmanharris",
        nperseg=4096,
        noverlap=4096 - 4096 // 4,
    )

    absed = np.abs(spectre)
    spectre = 20 * np.log10(absed / np.max(absed))
    spectre = np.where(spectre < -200, -200, spectre)

    plt.figure()
    plt.pcolormesh(time, freq, spectre, cmap="magma", shading="gouraud")
    plt.yscale("log")
    plt.ylim((20, samplerate / 2))
    plt.ylabel("Frequency [Hz]")
    plt.xlabel("Time [s]")
    plt.savefig("img/" + name + ".png")

Path("img").mkdir(parents=True, exist_ok=True)
for path in Path("snd").glob("chirp_*.wav"):
    print(f"Plotting {path}")
    data, samplerate = soundfile.read(str(path))
    plot(data, samplerate, path.stem)
    gc.collect()  # It goes out of memory with 32bit CPython (x86) on Windows.
