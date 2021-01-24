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

datapaths = [
    Path("snd/chirp_Simple.wav"),
    Path("snd/chirp_Mipmap.wav"),
    Path("snd/chirp_Lpsosc.wav"),
    Path("snd/chirp_Cpsosc.wav"),
]

Path("img").mkdir(parents=True, exist_ok=True)
for path in datapaths:
    data, samplerate = soundfile.read(str(path))
    plot(data, samplerate, path.stem)
# plt.show()
