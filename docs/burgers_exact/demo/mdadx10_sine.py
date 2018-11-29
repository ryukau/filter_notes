"""
Reference

- [NLALectureVT - NLALectureVT.pdf](http://perso.univ-lemans.fr/~vtournat/wa_files/NLALectureVT.pdf)
"""

import numpy
import matplotlib.pyplot as pyplot
import pathlib
import soundfile

from pyfftw.interfaces.numpy_fft import fft, rfft, hfft
from scipy.special import iv


def to_decibel(data):
    return 20 * numpy.log10(abs(data / numpy.max(abs(data))))


class Plotter:
    def __init__(self, out_dir, num_freq_bin):
        self.out_dir = out_dir
        self.index = 0

    def plot(self, data, freq, x, w):
        dpi = 96
        fig = pyplot.figure(figsize=(1280 / dpi, 720 / dpi), dpi=dpi)

        ax1 = fig.add_subplot(121)
        ax1.set_title("Time Domain")
        ax1.grid()
        ax1.set_xlabel("Phase[rad/π]")
        ax1.set_ylabel("Pressure")
        ax1.set_ylim((-1.1, 1.1))
        ax1.plot(x, data)

        ax2 = fig.add_subplot(122)
        ax2.set_title("Power Frequency")
        ax2.grid()
        ax2.set_xlabel("Harmonics")
        ax2.set_ylabel("Gain[dB]")
        ax2.set_ylim((-120, 10))
        ax2.plot(freq)

        fig.suptitle(
            "Normalized 5th Order Sine Approximation from mdaDX10.cpp , w={:f}\nx + x * x * x * (w * x * x - 1.0 - w)".
            format(w))
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)

        # pyplot.show()
        fig.savefig(
            self.out_dir + "/mdadx10_sine_{:04d}.png".format(self.index))
        pyplot.close()

        print(self.index)
        self.index += 1


def mdadx10_sine(x, w):
    return x + x * x * x * (w * x * x - 1.0 - w)


def normalize(data):
    max_value = numpy.max(numpy.abs(data))
    if max_value == 0:
        return data
    return data / max_value


def render_plot():
    out_dir = "img"
    pathlib.Path(out_dir).mkdir(exist_ok=True)

    x = numpy.linspace(-1, 1, 512)
    rich = numpy.linspace(0, 1, 1201)
    rich = 0.5 - 3 * rich * rich

    plotter = Plotter(out_dir, int(len(x) / 2) + 1)

    for w in rich:
        data = normalize(mdadx10_sine(x, w))
        freq = to_decibel(rfft(data, planner_effort='FFTW_ESTIMATE'))
        plotter.plot(data, freq, x, w)


def render_sound():
    """
    60 Hz で音を作る。 44100 / 60 = 735
    """
    x = numpy.linspace(-1, 1, 735)
    rich = numpy.linspace(0, 1, 1201)
    rich = 0.5 - 3 * rich * rich

    wave = []
    for w in rich:
        wave.append(normalize(mdadx10_sine(x, w)))

    data = numpy.concatenate(wave)
    soundfile.write("mdadx10_sine.wav", data, 44100, subtype="FLOAT")


if __name__ == "__main__":
    render_plot()
    render_sound()
