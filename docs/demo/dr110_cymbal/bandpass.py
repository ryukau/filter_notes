import matplotlib.pyplot as plt
import numpy
import soundfile
import scipy.signal as signal


def plot_frequency_response(systems):
    dpi = 96
    fig = plt.figure(figsize=(640 / dpi, 800 / dpi), dpi=dpi)

    ax1 = fig.add_subplot(211)
    ax1.set_title("Amplitude")
    ax1.set_xlabel("Frequency[Hz]")
    ax1.set_ylabel("Gain[dB]")
    ax1.grid(which="both")

    ax2 = fig.add_subplot(212)
    ax2.set_title("Phase")
    ax2.set_xlabel("Frequency[Hz]")
    ax2.set_ylabel("Phase[deg]")
    ax2.set_ylim((-200, 200))
    ax2.grid(which="both")

    freq = numpy.geomspace(10, 1e5, 256)
    for name, system in systems.items():
        w, mag, phase = signal.bode(system, freq)
        ax1.semilogx(w / 2 / numpy.pi, mag, label=name)  # Bode magnitude plot
        ax2.semilogx(w / 2 / numpy.pi, phase, label=name)  # Bode phase plot
    ax1.legend()
    ax2.legend()

    fig.suptitle("Bode Plot of DR-110 Band-pass Filter")
    fig.tight_layout()
    fig.subplots_adjust(top=0.92)

    plt.show()


def write_wav(samplerate, source, system, prefix=""):
    num, den, dt = signal.cont2discrete(system, 1.0 / samplerate)
    dest = signal.lfilter(num[0], den, source)

    soundfile.write(prefix + "_source.wav", source, samplerate)
    soundfile.write(prefix + "_dest.wav", dest, samplerate)


if __name__ == "__main__":
    samplerate = 44100

    source = numpy.random.uniform(-0.1, 0.1, samplerate * 2)

    sys_bp0 = (
        [9.471e+10, 2.5625e+16, 0],
        [3437973.0, 1.816815e+11, 8.03e+15, 3.125e+20],
    )

    sys_bp1 = (
        [2.68345e+11, 3.5234375e+16, 0],
        [30108309.0, 5.227075e+11, 1.56671875e+16, 1.953125e+20],
    )

    plot_frequency_response({"BP0": sys_bp0, "BP1": sys_bp1})
    # write_wav(samplerate, source, sys_bp0, "BP0")
    # write_wav(samplerate, source, sys_bp1, "BP1")
