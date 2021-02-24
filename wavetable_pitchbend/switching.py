import gc
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile
from pathlib import Path

def toDecibel(data):
    absed = np.abs(data)
    return 20 * np.log10(absed / np.max(absed))

def frequencyToMidinote(freq):
    return 12 * np.log2(freq / 440) + 69

def midinoteToFrequency(note):
    return 440 * np.exp2((note - 69) / 12)

def saw_spectrum(size):
    spec = [0] + [(-1)**k / k for k in range(1, int(size / 2 + 1))]
    spec = -1j * size / np.pi * np.array(spec)
    return spec

def tri_spectrum(size):
    spec = np.zeros(size / 2 + 1)
    for n in range(1, len(spec), 2):
        spec[n] = (-1)**((n - 1) / 2) / (n * n)
    return spec * 8 / np.pi**2

def plotSos(sos):
    w, h = signal.sosfreqz(sos)
    amp = toDecibel(np.abs(h))
    plt.plot(w, amp)
    plt.grid()
    plt.show()

class TableOsc:
    def __init__(self, samplerate, oversample=2):
        self.fs = oversample * samplerate
        self.phase = 0

        exponent = int(np.log2(self.fs / 10))
        size = 2**exponent  # Lowest to 10 Hz.
        spectrum = saw_spectrum(size)
        specSize = size / 2

        self.basefreq = self.fs / (2 * size)
        self.basenote = frequencyToMidinote(self.basefreq)

        self.table = []
        for idx in range(exponent + 1):
            spec = np.zeros_like(spectrum)
            cutoff = int(specSize / 2**(idx)) + 1  # +1 for DC component.
            spec[:cutoff] = spectrum[:cutoff]
            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl, tbl[0]))
            self.table.append(tbl)

        self.prevOct = 0  # debug

    def process(self, note):
        """
        note is midinote number.
        """
        self.phase += midinoteToFrequency(note) / self.fs
        self.phase -= np.floor(self.phase)

        octave = int(np.clip((note - self.basenote) / 12, 0, len(self.table) - 1))
        pos = (len(self.table[0]) - 1) * self.phase
        idx = int(pos)
        frac = pos - idx

        # debug
        isTableSwitched = 1 if self.prevOct != octave else 0
        self.prevOct = octave

        x0 = self.table[octave][idx]
        x1 = self.table[octave][idx + 1]
        return (x0 + frac * (x1 - x0), isTableSwitched)

class TableOscBilinear:
    def __init__(self, samplerate, oversample=2):
        self.fs = oversample * samplerate
        self.phase = 0

        exponent = int(np.log2(self.fs / 10))
        size = 2**exponent  # Lowest to 10 Hz.
        spectrum = saw_spectrum(size)
        specSize = size / 2

        self.basefreq = self.fs / (2 * size)
        self.basenote = frequencyToMidinote(self.basefreq)

        self.table = []
        for idx in range(exponent + 1):
            spec = np.zeros_like(spectrum)
            cutoff = int(specSize / 2**(idx)) + 1  # +1 for DC component.
            spec[:cutoff] = spectrum[:cutoff]
            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl, tbl[0]))
            self.table.append(tbl)

        self.prevOct = 0  # debug

    def process(self, note):
        """
        note is midinote number.
        """
        self.phase += midinoteToFrequency(note) / self.fs
        self.phase -= np.floor(self.phase)

        octFloat = np.clip((note - self.basenote) / 12, 0, len(self.table) - 2)
        octave = int(octFloat)
        yFrac = octFloat - octave

        pos = (len(self.table[0]) - 1) * self.phase
        idx = int(pos)
        xFrac = pos - idx

        # debug
        isTableSwitched = 1 if self.prevOct != octave else 0
        self.prevOct = octave

        x0 = self.table[octave][idx]
        x1 = self.table[octave][idx + 1]
        s0 = x0 + xFrac * (x1 - x0)

        x0 = self.table[octave + 1][idx]
        x1 = self.table[octave + 1][idx + 1]
        s1 = x0 + xFrac * (x1 - x0)

        return (s0 + yFrac * (x1 - s0), isTableSwitched)

def plotSpectrogram(sig, debugSig, samplerate, name):
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

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches(8, 6)
    fig.set_tight_layout(True)
    ax[0].pcolormesh(time, freq, spectre, cmap="magma", shading="gouraud")
    ax[0].set_yscale("log")
    ax[0].set_ylim((20, samplerate / 2))
    ax[0].set_ylabel("Frequency [Hz]")
    ax[0].set_xlabel("Time [s]")

    time = np.linspace(0, len(debugSig) / samplerate, len(debugSig))
    ax[1].plot(time, debugSig, lw=1, color="black")
    ax[1].grid()

    plt.savefig("img/" + name + ".png")
    plt.close("all")

def testOsc(Osc, oversample: int):
    print(f"Processing: {oversample}x {Osc.__name__}")

    gain = 0.25
    samplerate = 48000
    lowpass = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=oversample)
    # plotSos(lowpass)

    osc = Osc(samplerate, oversample)

    sig = np.empty(8 * oversample * samplerate)
    switched = np.empty_like(sig)  # debug

    note = np.linspace(0, 128, len(sig) // 2)
    note = np.hstack((note, np.flip(note)))
    for i in range(len(sig)):
        sig[i], switched[i] = osc.process(note[i])
    sig = gain * signal.sosfilt(lowpass, sig)[::oversample]
    switched = signal.sosfilt(lowpass, switched)[::oversample]  # debug

    Path("snd").mkdir(parents=True, exist_ok=True)

    name = f"{Osc.__name__}_x{oversample}_switch"
    soundfile.write(f"snd/{name}.wav", sig, samplerate, subtype="FLOAT")

    plotSpectrogram(sig, switched, samplerate, name)  # debug

    gc.collect()  # Maybe goes out of memory with 32bit CPython.

testOsc(TableOsc, 2)
testOsc(TableOscBilinear, 2)
testOsc(TableOsc, 4)
testOsc(TableOscBilinear, 4)
testOsc(TableOsc, 8)
testOsc(TableOscBilinear, 8)
