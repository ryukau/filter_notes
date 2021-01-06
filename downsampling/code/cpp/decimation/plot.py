import matplotlib.pyplot as plt
import numpy as np
import soundfile
from pathlib import Path
from multiprocessing import Pool

def getPowerSpectrum(data):
    spectrum = np.fft.rfft(data)
    absed = np.abs(spectrum)
    power = 20 * np.log10(absed / np.max(absed))
    return power

def job(args):
    path, additivePower = args
    data, samplerate = soundfile.read(str(path), always_2d=True)
    data = data.T[0]
    power = getPowerSpectrum(data)
    freq = np.fft.rfftfreq(len(data), 1 / samplerate)

    plt.figure(figsize=(12, 8), tight_layout=True)
    plt.title(path.stem)
    plt.plot(freq, power, color="black", lw=1, label=path.stem)
    plt.plot(freq, additivePower, color="red", lw=1, ls=":", label="Additive")
    plt.grid()
    plt.legend()
    plt.ylim((-100, 5))
    plt.savefig(f"img/{path.stem}.png", dpi=100)

if __name__ == "__main__":
    imgDir = Path("img")
    if not imgDir.exists():
        imgDir.mkdir(parents=True, exist_ok=True)

    data, samplerate = soundfile.read("sample/saw1.wav", always_2d=True)
    additivePower = getPowerSpectrum(data.T[0])

    with Pool() as pool:
        args = [(path, additivePower) for path in Path("snd").glob("*.wav")]
        for _ in pool.imap_unordered(job, args):
            pass
