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

        self.basefreq = self.fs / size
        self.basenote = frequencyToMidinote(self.basefreq)

        self.table = []
        for idx in range(2, exponent + 1):
            spec = np.zeros_like(spectrum)
            cutoff = int(size / 2**(idx))
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

        self.basefreq = self.fs / size
        self.basenote = frequencyToMidinote(self.basefreq)

        self.table = []
        for idx in range(2, exponent + 1):
            spec = np.zeros_like(spectrum)
            cutoff = int(size / 2**(idx))
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

        return s0 + yFrac * (x1 - s0)

class MipmapOsc:
    def __init__(self, samplerate, oversample=2):
        self.fs = oversample * samplerate
        self.phase = 0

        exponent = int(np.log2(self.fs / 10))  # Lowest to 10 Hz.
        size = 2**exponent
        spectrum = saw_spectrum(size)

        self.basefreq = self.fs / size
        self.basenote = frequencyToMidinote(self.basefreq)

        nOctave = exponent - 2
        self.table = []
        for idx in range(nOctave):
            cutoff = int(size / 2**(idx + 2)) + 1
            spec = spectrum[:cutoff] / 2**(idx + 1)
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

class LpsOsc:
    """
    Wavetable oscillator used in CubicPadSynth. Bilinear interpolation.
    """
    def __init__(self, samplerate):
        self.fs = samplerate
        self.phase = 0

        minFreq = midinoteToFrequency(0)
        expFloat = np.log2(self.fs / np.floor(minFreq))
        exponent = int(np.ceil(expFloat))
        self.size = 2**exponent
        spectrum = saw_spectrum(self.size)

        self.table = []
        self.maxNote = np.ceil(frequencyToMidinote(20000))
        spec = np.zeros_like(spectrum)
        for idx in range(int(self.maxNote)):
            freq = midinoteToFrequency(idx)
            cutoff = int(len(spectrum) * minFreq / freq)

            spec[:cutoff] = spectrum[:cutoff]
            spec[cutoff:] = 0

            tbl = np.fft.irfft(spec)
            tbl = np.hstack((tbl, tbl[0]))
            self.table.append(tbl)
        self.table.append(np.zeros_like(self.table[0]))

    def processPhase(self, note):
        tick = midinoteToFrequency(note) * self.size / self.fs
        if tick >= self.size or tick < 0:
            tick = 0

        self.phase += tick
        if self.phase >= self.size:
            self.phase -= self.size

    def processLow(self, note):
        idx = int(self.phase)
        a0 = self.table[0][idx]
        a1 = self.table[0][idx + 1]
        fracX = self.phase - np.floor(self.phase)
        return a0 + fracX * (a1 - a0)

    def process(self, note):
        """
        note is midinote number.
        """
        self.processPhase(note)

        if note < 0:
            return processLow(note)

        nn = int(note)
        if nn >= self.maxNote:
            nn = self.maxNote - 1

        idx = int(self.phase)
        a0 = self.table[nn][idx]
        a1 = self.table[nn][idx + 1]
        b0 = self.table[nn + 1][idx]
        b1 = self.table[nn + 1][idx + 1]

        fracX = self.phase - np.floor(self.phase)
        x0 = a0 + fracX * (a1 - a0)
        x1 = b0 + fracX * (b1 - b0)

        fracY = note - nn
        return x0 + fracY * (x1 - x0)

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

class CpsOsc:
    """
    Wavetable oscillator used in CubicPadSynth.
    """
    def __init__(self, samplerate):
        self.fs = samplerate
        self.phase = 0

        minFreq = midinoteToFrequency(0)
        expFloat = np.log2(self.fs / np.floor(minFreq))
        exponent = int(np.ceil(expFloat))
        self.size = 2**exponent
        spectrum = saw_spectrum(self.size)

        self.table = []
        self.maxNote = np.ceil(frequencyToMidinote(20000))
        spec = np.zeros_like(spectrum)
        for idx in range(int(self.maxNote)):
            freq = midinoteToFrequency(idx)
            cutoff = int(len(spectrum) * minFreq / freq)

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

def plotSpectrum(sig):
    spec = toDecibel(np.abs(np.fft.rfft(sig)))
    plt.plot(spec, label="True", alpha=0.75, lw=2, color="blue")
    plt.grid()
    plt.legend()
    plt.ylim((-100, 10))
    plt.show()

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

def testTableOsc():
    gain = 0.25
    samplerate = 48000
    lowpass = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=2)
    # plotSos(lowpass)

    osc = TableOsc(samplerate)

    sig = np.empty(16 * samplerate)
    note = np.linspace(0, 128, len(sig) // 2)
    note = np.hstack((note, np.flip(note)))
    for i in range(len(sig)):
        sig[i] = osc.process(note[i])
    sig = gain * signal.sosfilt(lowpass, sig)[::2]

    Path("snd").mkdir(parents=True, exist_ok=True)
    soundfile.write("snd/simple.wav", sig, samplerate, subtype="FLOAT")

    # plotSpectrum(sig)
    plotSpectrogram(sig, samplerate, "TableOsc")

def testTableOscBilinear():
    gain = 0.25
    samplerate = 48000
    lowpass = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=2)
    # plotSos(lowpass)

    osc = TableOscBilinear(samplerate)

    sig = np.empty(16 * samplerate)
    note = np.linspace(0, 128, len(sig) // 2)
    note = np.hstack((note, np.flip(note)))
    for i in range(len(sig)):
        sig[i] = osc.process(note[i])
    sig = gain * signal.sosfilt(lowpass, sig)[::2]

    Path("snd").mkdir(parents=True, exist_ok=True)
    soundfile.write("snd/bilinear.wav", sig, samplerate, subtype="FLOAT")

    # plotSpectrum(sig)
    plotSpectrogram(sig, samplerate, "TableOscBilinear")

def testMipmapOsc():
    gain = 0.25
    samplerate = 48000
    lowpass = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=2)

    osc = MipmapOsc(samplerate)

    sig = np.empty(16 * samplerate)
    note = np.linspace(0, 128, len(sig) // 2)
    note = np.hstack((note, np.flip(note)))
    for i in range(len(sig)):
        sig[i] = osc.process(note[i])
    sig = gain * signal.sosfilt(lowpass, sig)[::2]

    Path("snd").mkdir(parents=True, exist_ok=True)
    soundfile.write("snd/mipmap.wav", sig, samplerate, subtype="FLOAT")

    # plotSpectrum(sig)
    plotSpectrogram(sig, samplerate, "MipmapOsc")

def testLpsOsc():
    gain = 0.25
    samplerate = 48000

    osc = LpsOsc(samplerate)

    sig = np.empty(8 * samplerate)
    note = np.linspace(0, 128, len(sig) // 2)
    note = np.hstack((note, np.flip(note)))
    for i in range(len(sig)):
        sig[i] = osc.process(note[i])
    sig = gain * sig

    Path("snd").mkdir(parents=True, exist_ok=True)
    soundfile.write("snd/lpsosc.wav", sig, samplerate, subtype="FLOAT")

    # plotSpectrum(sig)
    plotSpectrogram(sig, samplerate, "LpsOsc")

def testCpsOsc():
    gain = 0.25
    samplerate = 48000

    osc = CpsOsc(samplerate)

    sig = np.empty(8 * samplerate)
    note = np.linspace(0, 128, len(sig) // 2)
    note = np.hstack((note, np.flip(note)))
    for i in range(len(sig)):
        sig[i] = osc.process(note[i])
    sig = gain * sig

    Path("snd").mkdir(parents=True, exist_ok=True)
    soundfile.write("snd/cpsosc.wav", sig, samplerate, subtype="FLOAT")

    # plotSpectrum(sig)
    plotSpectrogram(sig, samplerate, "CpsOsc")

testTableOsc()
testTableOscBilinear()
testMipmapOsc()
testLpsOsc()
testCpsOsc()
