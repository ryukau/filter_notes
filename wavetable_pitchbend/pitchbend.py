import gc
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import soundfile
from pathlib import Path

def frequencyToMidinote(freq):
    return 12 * np.log2(freq / 440) + 69

def midinoteToFrequency(note):
    return 440 * np.exp2((note - 69) / 12)

def saw_spectrum(size):
    spec = [0] + [(-1)**k / k for k in range(1, int(size / 2 + 1))]
    spec = -1j * size / np.pi * np.array(spec)
    return spec

def saw_spectrum_shifted(size):
    spec = [0] + [-1 / k for k in range(1, int(size / 2 + 1))]
    spec = (-1j * size / np.pi) * np.array(spec)
    return spec

def tri_spectrum(size):
    spec = np.zeros(int(size / 2 + 1))
    for n in range(1, len(spec), 2):
        sign = -1 if (n - 1) % 2 == 0 else 1
        spec[n] = sign / (n * n)
    return spec * (size * 4 / np.pi**2)

def cubicInterp(y0, y1, y2, y3, t):
    """
    Range of t is in [0, 1]. Interpolates in between y1 and y2.
    """
    t2 = t * t
    c0 = y1 - y2
    c1 = (y2 - y0) / 2
    c2 = c0 + c1
    c3 = c0 + c2 + (y3 - y1) / 2
    return c3 * t * t2 - (c2 + c3) * t2 + c1 * t + y1

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

        x0 = self.table[octave][idx]
        x1 = self.table[octave][idx + 1]
        return x0 + frac * (x1 - x0)

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

        x0 = self.table[octave][idx]
        x1 = self.table[octave][idx + 1]
        s0 = x0 + xFrac * (x1 - x0)

        x0 = self.table[octave + 1][idx]
        x1 = self.table[octave + 1][idx + 1]
        s1 = x0 + xFrac * (x1 - x0)

        return s0 + yFrac * (s1 - s0)

class TableOscAltInterval:
    def __init__(self, samplerate, oversample=2):
        self.fs = oversample * samplerate
        self.phase = 0

        size = 2**int(np.log2(self.fs / 10))  # Lowest to 10 Hz.
        specSize = size / 2
        spectrum = saw_spectrum(size)

        bendRange = np.sqrt(3)
        nTable = int(-np.log(1 / specSize) / np.log(bendRange))

        self.basefreq = self.fs / size
        self.basenote = frequencyToMidinote(self.basefreq)
        self.interval = 12 * np.log2(bendRange)

        self.table = []
        for idx in range(nTable):
            spec = np.zeros_like(spectrum)
            cutoff = int(specSize * bendRange**(-idx))
            spec[:cutoff] = spectrum[:cutoff]
            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl, tbl[0]))
            self.table.append(tbl)
        self.table.append(np.zeros_like(self.table[0]))

    def process(self, note):
        """
        note is midinote number.
        """
        self.phase += midinoteToFrequency(note) / self.fs
        self.phase -= np.floor(self.phase)

        octFloat = np.clip((note - self.basenote) / self.interval, 0, len(self.table) - 2)
        iTbl = int(octFloat)
        yFrac = octFloat - iTbl

        pos = (len(self.table[0]) - 1) * self.phase
        idx = int(pos)
        xFrac = pos - idx

        x0 = self.table[iTbl][idx]
        x1 = self.table[iTbl][idx + 1]
        s0 = x0 + xFrac * (x1 - x0)

        x0 = self.table[iTbl + 1][idx]
        x1 = self.table[iTbl + 1][idx + 1]
        s1 = x0 + xFrac * (x1 - x0)

        return s0 + yFrac * (s1 - s0)

class MipmapOsc:
    def __init__(self, samplerate, oversample=2):
        self.fs = oversample * samplerate
        self.phase = 0

        exponent = int(np.log2(self.fs / 10))  # Lowest to 10 Hz.
        size = 2**exponent
        spectrum = saw_spectrum(size)

        self.basefreq = self.fs / size
        self.basenote = frequencyToMidinote(self.basefreq)

        nOctave = exponent
        self.table = []
        for idx in range(nOctave):
            cutoff = int(size / 2**(idx + 1)) + 1
            spec = spectrum[:cutoff] / 2**idx
            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl, tbl[0]))
            self.table.append(tbl)

    def process(self, note):
        """
        note is midinote number.
        """
        self.phase += midinoteToFrequency(note) / self.fs
        self.phase -= np.floor(self.phase)

        octave = int(np.clip((note - self.basenote) / 12, 0, len(self.table) - 1))
        pos = (len(self.table[octave]) - 1) * self.phase
        idx = int(pos)
        frac = pos - idx

        x0 = self.table[octave][idx]
        x1 = self.table[octave][idx + 1]
        return x0 + frac * (x1 - x0)

class MipmapOscCubic:
    def __init__(self, samplerate, oversample=2):
        self.fs = oversample * samplerate
        self.phase = 0

        exponent = int(np.log2(self.fs / 10))  # Lowest to 10 Hz.
        size = 2**exponent
        spectrum = saw_spectrum(size)

        self.basefreq = self.fs / size
        self.basenote = frequencyToMidinote(self.basefreq)

        nOctave = exponent
        self.table = []
        for idx in range(nOctave):
            cutoff = int(size / 2**(idx + 1)) + 1
            spec = spectrum[:cutoff] / 2**idx
            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl[-1], tbl, tbl[0], tbl[1]))
            self.table.append(tbl)
        self.ylim = len(self.table) - 1

    def process(self, note):
        """
        note is midinote number.
        """
        self.phase += midinoteToFrequency(note) / self.fs
        self.phase -= np.floor(self.phase)

        octave = int(np.clip((note - self.basenote) / 12, 0, self.ylim))
        pos = (len(self.table[octave]) - 3) * self.phase
        idx = int(pos)
        frac = pos - idx

        return cubicInterp(
            self.table[octave][idx],
            self.table[octave][idx + 1],
            self.table[octave][idx + 2],
            self.table[octave][idx + 3],
            pos - idx,
        )

class LpsOsc:
    """
    Wavetable oscillator used in CubicPadSynth. Bilinear interpolation.
    """
    def __init__(self, samplerate):
        self.fs = samplerate
        self.phase = 0

        expFloat = np.log2(self.fs / np.floor(midinoteToFrequency(0)))
        exponent = int(np.ceil(expFloat))
        self.size = 2**exponent
        minFreq = self.fs / self.size
        spectrum = saw_spectrum(self.size)

        self.table = []
        self.maxNote = np.ceil(frequencyToMidinote(20000))
        spec = np.zeros_like(spectrum)
        for idx in range(int(self.maxNote)):
            freq = midinoteToFrequency(idx)
            cutoff = int(self.size / 2 * minFreq / freq) + 1

            spec[:cutoff] = spectrum[:cutoff]
            spec[cutoff:] = 0

            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl, tbl[0]))
            self.table.append(tbl)
        self.table.append(np.zeros_like(self.table[0]))

    def process(self, note):
        """
        note is midinote number.
        """
        self.phase += midinoteToFrequency(note) / self.fs
        self.phase -= np.floor(self.phase)

        note = np.clip(note, 0, self.size - 2)
        nn = int(note)

        phs = self.size * self.phase
        idx = int(phs)
        a0 = self.table[nn][idx]
        a1 = self.table[nn][idx + 1]
        b0 = self.table[nn + 1][idx]
        b1 = self.table[nn + 1][idx + 1]

        fracX = phs - idx
        x0 = a0 + fracX * (a1 - a0)
        x1 = b0 + fracX * (b1 - b0)

        fracY = note - nn
        return x0 + fracY * (x1 - x0)

class CpsOsc:
    """
    Wavetable oscillator used in CubicPadSynth.
    """
    def __init__(self, samplerate):
        self.fs = samplerate
        self.phase = 0

        expFloat = np.log2(self.fs / np.floor(midinoteToFrequency(0)))
        exponent = int(np.ceil(expFloat))
        self.size = 2**exponent
        minFreq = self.fs / self.size
        spectrum = saw_spectrum(self.size)

        self.table = []
        self.maxNote = np.ceil(frequencyToMidinote(20000))
        spec = np.zeros_like(spectrum)
        for idx in range(int(self.maxNote)):
            freq = midinoteToFrequency(idx)
            cutoff = int((len(spectrum) - 1) * minFreq / freq) + 1

            spec[:cutoff] = spectrum[:cutoff]
            spec[cutoff:] = 0

            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl[-1], tbl, tbl[0], tbl[1]))
            self.table.append(tbl)
        self.table.append(np.zeros_like(self.table[0]))

    def processPhase(self, note):
        tick = midinoteToFrequency(note) * self.size / self.fs
        if tick >= self.size or tick < 0:
            tick = 0

        self.phase += tick
        if self.phase >= self.size:
            self.phase -= self.size

    def process(self, note):
        """
        note is midinote number.
        """
        self.processPhase(note)

        nn = int(note)
        if nn < 0:
            nn = 0
        elif nn > self.maxNote - 1:
            nn = self.maxNote - 1

        idx = int(self.phase)
        return cubicInterp(
            self.table[nn][idx],
            self.table[nn][idx + 1],
            self.table[nn][idx + 2],
            self.table[nn][idx + 3],
            self.phase - np.floor(self.phase),
        )

def plotSpectrogram(sig, samplerate, name):
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

    plt.figure(figsize=(8, 4), tight_layout=True)
    plt.pcolormesh(time, freq, spectre, cmap="magma", shading="gouraud")
    plt.yscale("log")
    plt.ylim((20, samplerate / 2))
    plt.ylabel("Frequency [Hz]")
    plt.xlabel("Time [s]")
    plt.savefig("img/" + name + ".png")
    plt.close("all")

    gc.collect()  # Maybe out of memory on 32bit CPython.

def testOsc(Osc, oversample, duration=8):
    print(f"Processing {Osc.__name__}")

    gain = 0.25
    samplerate = 48000

    osc = Osc(samplerate)

    sig = np.empty(oversample * duration * samplerate)
    note = np.linspace(0, 128, len(sig) // 2)
    note = np.hstack((note, np.flip(note)))
    for i in range(len(sig)):
        sig[i] = osc.process(note[i])

    if oversample >= 2:
        lowpass = signal.ellip(12, 0.01, 100, 0.48, "low", output="sos", fs=oversample)
        sig = signal.sosfilt(lowpass, sig)[::oversample]

    sig *= gain

    Path("snd").mkdir(parents=True, exist_ok=True)
    soundfile.write(f"snd/chirp_{Osc.__name__}.wav", sig, samplerate, subtype="FLOAT")

    plotSpectrogram(sig, samplerate, Osc.__name__)

testOsc(LpsOsc, 1)
testOsc(CpsOsc, 1)
testOsc(TableOsc, 2)
testOsc(TableOscBilinear, 2)
testOsc(TableOscAltInterval, 2)
testOsc(MipmapOsc, 2)
testOsc(MipmapOscCubic, 2)
