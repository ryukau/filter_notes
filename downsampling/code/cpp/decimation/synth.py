import numpy as np
import soundfile
from pathlib import Path
from multiprocessing import Pool

def generateSaw(samplerate, duration, frequency):
    time = np.linspace(0, duration, int(samplerate * duration))
    phase = frequency * time

    omega_t = 2 * np.pi * phase
    sig = np.zeros_like(omega_t)
    overtone = np.arange(frequency, samplerate / 2, frequency)
    for k, freq in enumerate(overtone, 1):
        sig += ((-1)**k / k) * np.sin(k * omega_t)
    return 2 * sig / np.pi

def writeSaw(outDir, oversample, samplerate, duration, frequency):
    fs = oversample * samplerate
    data = generateSaw(fs, 1, frequency)
    path = outDir / f"saw{oversample}.wav"
    soundfile.write(str(path), data, fs, subtype="FLOAT")

def midinoteToFreq(semitone):
    return 440 * np.power(2, (semitone - 69) / 12)

if __name__ == "__main__":
    sampleDir = Path("sample")
    if not sampleDir.exists():
        sampleDir.mkdir(parents=True, exist_ok=True)

    samplerate = 48000
    note = 84
    writeSaw(sampleDir, 1, samplerate, 1, midinoteToFreq(note))
    writeSaw(sampleDir, 2, samplerate, 1, midinoteToFreq(note))
    writeSaw(sampleDir, 4, samplerate, 1, midinoteToFreq(note))
    writeSaw(sampleDir, 8, samplerate, 1, midinoteToFreq(note))
    writeSaw(sampleDir, 16, samplerate, 1, midinoteToFreq(note))
