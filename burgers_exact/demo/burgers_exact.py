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


def burgers_hopf_cole(xi, theta, gamma, loop=128):
    """
    p_a / p_0 を返す。
    """
    gamma_div_2 = gamma / 2
    xi_div_gamma = xi / gamma

    numer = 0
    denom = 0
    sign = 1
    for n in range(1, loop + 1):
        iv_exp = iv(n, gamma_div_2) * numpy.exp(-n * n * xi_div_gamma)
        n_theta = n * theta
        numer += sign * n * iv_exp * numpy.sin(n_theta)
        sign = -sign
        denom += sign * iv_exp * numpy.cos(n_theta)
    return (4 / gamma * numer) / (iv(0, gamma_div_2) + 2 * denom)


def to_decibel(data):
    return 20 * numpy.log10(abs(data / numpy.max(abs(data))))


class Plotter:
    def __init__(self, out_dir, num_freq_bin):
        self.out_dir = out_dir
        self.index = 0

        self.time_xtick_major = numpy.linspace(0, 2 * numpy.pi, 9)
        self.time_ytick_major = numpy.linspace(-1, 1, 9)

    def plot(self, data, freq, theta, xi, gamma):
        dpi = 96
        fig = pyplot.figure(figsize=(1280 / dpi, 720 / dpi), dpi=dpi)

        ax1 = fig.add_subplot(121)
        ax1.set_title("Time Domain")
        ax1.grid()
        ax1.set_xticks(self.time_xtick_major)
        ax1.set_yticks(self.time_ytick_major)
        ax1.set_xlabel("Phase[rad]")
        ax1.set_ylabel("Pressure")
        ax1.set_ylim((-1.1, 1.1))
        ax1.plot(theta, data)

        ax2 = fig.add_subplot(122)
        ax2.set_title("Power Frequency")
        ax2.grid()
        ax2.set_xlabel("Harmonics")
        ax2.set_ylabel("Gain[dB]")
        ax2.set_ylim((-120, 10))
        ax2.plot(freq)

        fig.suptitle("Solution of Burgers Equation, ξ={:f}, Γ={:d}".format(
            xi, gamma))
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)

        # pyplot.show()
        fig.savefig(self.out_dir + "/burgers_{:04d}.png".format(self.index))
        pyplot.close()

        print(self.index)
        self.index += 1


def render_plot():
    out_dir = "img"
    pathlib.Path(out_dir).mkdir(exist_ok=True)

    theta = numpy.linspace(0, 2 * numpy.pi, 512)
    gamma = numpy.arange(1, 31)

    plotter = Plotter(out_dir, int(len(theta) / 2) + 1)

    for gam in gamma:
        xi = numpy.linspace(0, 10, 121)
        for x in xi:
            data = burgers_hopf_cole(x, theta, gam, 128)
            freq = to_decibel(rfft(data, planner_effort='FFTW_ESTIMATE'))
            plotter.plot(data, freq, theta, x, gam)


def render_sound():
    """
    60 Hz で音を作る。 44100 / 60 = 735
    """
    theta = numpy.linspace(0, 2 * numpy.pi, 735, endpoint=False)
    gamma = numpy.arange(1, 31)

    wave = []
    for gam in gamma:
        xi = numpy.linspace(0, 10, 121)
        print(gam)
        for x in xi:
            wave.append(burgers_hopf_cole(x, theta, gam, 128))

    data = numpy.concatenate(wave)
    soundfile.write("burgers_exact.wav", data, 44100, subtype="FLOAT")


if __name__ == "__main__":
    render_plot()
    render_sound()
