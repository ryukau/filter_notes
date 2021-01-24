import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

def toDecibel(data):
    absed = np.abs(data)
    return 20 * np.log10(absed / np.max(absed))

def frequencyToMidinote(freq):
    return 12 * np.log2(freq / 440) + 69

def midinoteToFrequency(note):
    return 440 * np.exp2((note - 69) / 12)

def saw_spectrum_fft(size):
    wave = np.linspace(-1, 1, size)
    wave = np.roll(wave, size // 2)
    spec = np.fft.rfft(wave)
    return spec

def saw_spectrum(size):
    spec = [0] + [(-1)**k / k for k in range(1, int(size / 2 + 1))]
    spec = -1j * size / np.pi * np.array(spec)
    return spec

def plotTable(table, spectrum, fullSpectrum=None, fs=2, name="test"):
    fig, ax = plt.subplots(1, 2)
    fig.set_size_inches(8, 3)
    fig.set_tight_layout(True)

    ax[0].plot(table, color="black")
    ax[0].set_title("Wavetable")
    ax[0].set_xlabel("Time [sample]")
    ax[0].set_ylabel("Amplitude")

    freq = np.linspace(0, fs / 2, len(spectrum))
    spec = toDecibel(spectrum)
    spec = np.where(spec < -200, -200, spec)
    ax[1].plot(freq, spec, color="black", label="Bandlimited")
    ax[1].set_title("Spectrum")
    ax[1].set_xlabel(f"Frequency (fs={fs})")
    ax[1].set_xscale("log")
    ax[1].set_xlim((10, fs / 2))
    ax[1].set_ylabel("Amplitude [dB]")
    ax[1].set_ylim((-100, 5))

    if fullSpectrum is not None:
        spec = toDecibel(fullSpectrum)
        spec = np.where(spec < -200, -200, spec)
        ax[1].plot(freq, spec, alpha=0.5, color="red", ls="--", label="Source")
        ax[1].legend()

    for axis in ax:
        axis.grid()
    plt.savefig("img/" + name + ".svg")

def plotSpectrogram(sig, samplerate, name="test"):
    freq, time, spectre = signal.spectrogram(
        sig,
        samplerate,
        # ("tukey", 0.33),
        "blackmanharris",
        nperseg=512,
        noverlap=512 - 512 // 4,
    )

    absed = np.abs(spectre)
    spectre = 20 * np.log10(absed / np.max(absed))
    spectre = np.where(spectre < -100, -100, spectre)

    plt.figure(figsize=(4, 3.5), tight_layout=True)
    plt.title(name)
    plt.pcolormesh(time, freq, spectre, cmap="magma", shading="gouraud")
    # plt.yscale("log")
    # plt.ylim((20, samplerate / 2))
    plt.ylabel("Frequency [Hz]")
    plt.xlabel("Time [sec]")
    filename = name.replace(" ", "")
    plt.savefig("img/" + filename + ".png")

def renderTable(table, freqStart, freqEnd, fs):
    phase = np.linspace(freqStart / fs, freqEnd / fs, fs).cumsum() % 1.0
    x = np.linspace(0, 1, len(table))
    return np.interp(phase, x, table, period=1)

def testNaiveTable(fs=48000, f0=880, size=8192):
    source = saw_spectrum(size)

    spectrum = source.copy()
    # cutoff = int(fs / (2 * f0))
    # spectrum[cutoff:] = 0
    table = np.fft.irfft(spectrum)
    plotTable(table, spectrum, source, fs, "table")

    sigF0 = renderTable(table, f0, f0, fs)
    plotSpectrogram(sigF0, fs, "Constant Pitch")

    sigUp = renderTable(table, f0, 2 * f0, fs)
    plotSpectrogram(sigUp, fs, "Rising Pitch")

    sigDown = renderTable(table, f0, 0.5 * f0, fs)
    plotSpectrogram(sigDown, fs, "Falling Pitch")

    # plt.show()

def test2xOversampledTable(fs=48000, f0=880, size=8192):
    oversample = 2
    upRate = oversample * fs
    source = saw_spectrum(size)

    spectrum = source.copy()
    cutoff = int(upRate / (2 * oversample * f0))
    spectrum[cutoff:] = 0
    table = np.fft.irfft(spectrum)
    plotTable(table, spectrum, source, upRate, "table2x")

    sigUp = renderTable(table, f0, 2 * f0, upRate)
    plotSpectrogram(sigUp, upRate, "Rising Pitch 2x")

    lowpass = signal.ellip(12, 0.01, 100, 0.4, "low", output="sos", fs=2)
    sigFiltered = signal.sosfilt(lowpass, sigUp)
    plotSpectrogram(sigFiltered, upRate, "Filtered Rising Pitch 2x")

    sigDecimated = sigFiltered[::2]
    plotSpectrogram(sigDecimated, fs, "Filtered Rising Pitch 2x Decimated")

    # plt.show()

def plotTable(fs=48000, size=1024):
    oversample = 1
    source = saw_spectrum(size)

    nTable = 8

    spectrum = source.copy()
    table = []
    for idx in reversed(range(nTable)):
        f0 = fs // (2**(idx + 2))
        cutoff = int(fs / (2 * oversample * f0))
        spectrum[cutoff:] = 0
        table.append(np.fft.irfft(spectrum))
    table.reverse()

    fig, ax = plt.subplots(nTable, 1)
    fig.set_size_inches(4, 4)
    fig.set_tight_layout(True)
    for axis, tbl in zip(ax, table):
        axis.plot(tbl)
        axis.set_axis_off()
    plt.show()

# testNaiveTable()
# test2xOversampledTable()
plotTable()
