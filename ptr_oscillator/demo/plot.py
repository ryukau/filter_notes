import matplotlib.pyplot as pyplot
import numpy
import soundfile
import numpy.fft as fft

from scipy import signal
from pathlib import Path

def normalize(sig):
    max_value = numpy.max(numpy.abs(sig))
    if max_value == 0:
        return sig
    return sig / max_value

def to_decibel(data):
    data_abs = numpy.abs(data)
    return 20 * numpy.log10(data_abs / numpy.max(data_abs))

def plotSpectrum(waveform):
    nrow = 4
    ncol = 3
    fig, ax = pyplot.subplots(nrow, ncol)
    fig.suptitle(f"PTR{waveform}")
    fig.set_size_inches(12.8, 8.0)
    fig.tight_layout(pad=2, rect=[0, 0, 1, 0.95])

    ptr_wav = sorted(Path("snd").glob(f"PTR{waveform}[0-9]*.wav"))
    for index, wav in enumerate(ptr_wav):
        data, samplerate = soundfile.read(wav, always_2d=True)
        spec = to_decibel(numpy.abs(fft.rfft(data.T[0])))

        row = index % nrow
        col = index // nrow
        ax[row][col].set_title(f"W = {index}")
        ax[row][col].plot(spec, lw=1, color="black")
        ax[row][col].set_ylim((-100, 0))
        ax[row][col].grid()
    ax[nrow - 1][ncol - 1].axis("off")
    pyplot.savefig(f"img/{waveform.lower()}_spectrum.png", dpi=100)
    pyplot.clf()

def plotWave(waveform):
    cmap = pyplot.get_cmap("plasma")
    pyplot.figure(figsize=(9.6, 4.8))

    ptr_wav = sorted(Path("snd").glob(f"PTR{waveform}[0-9]*.wav"))
    for index, wav in enumerate(ptr_wav):
        if index < -1:
            continue
        data, samplerate = soundfile.read(wav, always_2d=True)
        sig = data.T[0][:1024]
        pyplot.plot(
            sig,
            label=wav.stem,
            alpha=0.75,
            lw=1,
            color=cmap(index / len(ptr_wav)),
        )
    pyplot.grid()
    pyplot.legend(loc=1)
    pyplot.xlim((0, 80))
    pyplot.ylim((-1.1, 1.1))
    pyplot.title(f"PTR{waveform} Waveform")
    pyplot.xlabel("Time [sample]")
    pyplot.ylabel("Amplitude")
    pyplot.savefig(f"img/{waveform.lower()}_waveform.png", dpi=100)
    pyplot.clf()

waveform = ["Saw", "Tri", "Ramp", "Step"]

for wf in waveform:
    plotSpectrum(wf)
    plotWave(wf)
