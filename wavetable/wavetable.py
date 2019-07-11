import itertools
import matplotlib.pyplot as pyplot
import numpy
import pprint
import scipy.signal as signal
import scipy.interpolate as interpolate
import soundfile

import numpy.fft as fft

def normalize(sig):
    max_value = numpy.max(numpy.abs(sig))
    if max_value == 0:
        return sig
    return sig / max_value

def to_decibel(data):
    data_abs = numpy.abs(data)
    return 20 * numpy.log10(data_abs / numpy.max(data_abs))

def compare(true_sig, real_sig):
    true_spec = to_decibel(numpy.abs(fft.rfft(true_sig)))
    real_spec = to_decibel(numpy.abs(fft.rfft(real_sig)))
    pyplot.plot(true_spec, label="True", alpha=0.75, lw=2, color="blue")
    pyplot.plot(real_spec, label="Real", alpha=0.75, lw=1, color="red")
    pyplot.grid()
    pyplot.legend()
    pyplot.ylim((-200, 10))
    pyplot.show()

def error_cubic_vs_sinc(true_sig, cubic, sinc):
    true_spec = to_decibel(numpy.abs(fft.rfft(true_sig)))
    spec1 = to_decibel(numpy.abs(fft.rfft(cubic))) - true_spec
    spec2 = to_decibel(numpy.abs(fft.rfft(sinc))) - true_spec
    pyplot.plot(spec1, label="Cubic", alpha=0.75, lw=1, color="red")
    pyplot.plot(spec2, label="Sinc", alpha=0.75, lw=1, color="green")
    pyplot.grid()
    pyplot.legend()
    pyplot.show()

def mean_absolute_error(true_sig, real_sig):
    true_spec = to_decibel(numpy.abs(fft.rfft(true_sig)))
    real_spec = to_decibel(numpy.abs(fft.rfft(real_sig)))
    return numpy.nanmean(numpy.abs(true_spec - real_spec))

def calc_error(additive, naive, linterp, cubic, sinc):
    return {
        "naive": mean_absolute_error(additive, naive),
        "linterp": mean_absolute_error(additive, linterp),
        "cubic": mean_absolute_error(additive, cubic),
        "sinc": mean_absolute_error(additive, sinc),
    }

def downsampling(samplerate, sig, ratio):
    nyquist = samplerate / 2
    cutoff = 21000 / nyquist / ratio if nyquist > 21000 else nyquist
    lowpass = signal.firwin(512, cutoff, window=("kaiser", 8))
    for _ in range(2):
        sig = signal.convolve(sig, lowpass, mode="same")
    return sig[::ratio]

def additive_saw(samplerate, base_freq, phase):
    """
    amp ~= [-1, 1]. はみ出す。
    """
    omega_t = 2 * numpy.pi * phase
    sig = numpy.zeros_like(omega_t)
    overtone = numpy.arange(base_freq, samplerate / 2, base_freq)
    for k, freq in enumerate(overtone, 1):
        sig += ((-1)**k / k) * numpy.sin(k * omega_t)
    return 2 * sig / numpy.pi

def saw_spectrum(size):
    spec = [0] + [(-1)**k / k for k in range(1, int(size / 2 + 1))]
    spec = -1j * size / numpy.pi * numpy.array(spec)
    return spec

def naive_saw(samplerate, phase, size=512):
    spec = saw_spectrum(size)

    table = numpy.fft.irfft(spec)
    table = numpy.append(table, [table[0]])

    xp = numpy.linspace(0, 1, len(table))
    return numpy.interp(phase % 1.0, xp, table)

def yoshimi_saw_linterp(samplerate, base_freq, phase, size=512):
    spec = saw_spectrum(size)
    n_harmonics = int(samplerate / 2 / base_freq)
    spec[n_harmonics + 1:] = 0

    table = numpy.fft.irfft(spec)
    table = numpy.append(table, [table[0]])

    xp = numpy.linspace(0, 1, len(table))
    return numpy.interp(phase % 1.0, xp, table)

def yoshimi_saw_cubic(samplerate, base_freq, phase, size=512):
    spec = saw_spectrum(size)
    n_harmonics = int(samplerate / 2 / base_freq)
    spec[n_harmonics + 1:] = 0

    table = numpy.fft.irfft(spec)
    table = numpy.append(table, [table[0]])

    xp = numpy.linspace(0, 1, len(table))
    interp = interpolate.CubicSpline(xp, table, bc_type="periodic")
    return interp(phase % 1.0)

def sinc_saw(samplerate, base_freq, phase, size=512):
    spec = saw_spectrum(size)
    n_harmonics = int(samplerate / 2 / base_freq)
    spec[n_harmonics + 1:] = 0

    table = numpy.fft.irfft(spec)

    window = signal.windows.kaiser(len(table) - 1, 2.0952)
    # window = signal.windows.kaiser(len(table) // 2 - 1, 2.0952)
    window = numpy.insert(window, 0, 0)
    half = len(window) // 2
    w_index = numpy.arange(-half, half)

    phase = phase % 1.0 * len(table)
    sig = numpy.zeros_like(phase)
    for sig_index, ph in enumerate(phase):
        # print(f"\r{sig_index} / {len(sig)}", end="")
        win = window * numpy.sinc(ph % 1.0 - w_index)
        rolled = numpy.roll(table, half - int(ph))[:len(window)]
        sig[sig_index] = numpy.sum(win * rolled)
    # print()
    return sig

def generate_wave(samplerate, duration, frequency, type="sin"):
    time = numpy.linspace(0, duration, int(samplerate * duration))
    phase = 2 * numpy.pi * frequency * time
    if type == "saw":
        return signal.sawtooth(phase, 1)
    if type == "square":
        return signal.square(phase, 0.5)
    return numpy.sin(phase)

def fm_phase(samplerate, base_freq, mod_amount, mod):
    """
    freq = base ± lfo_amount.
    """
    freq = base_freq + mod_amount * mod
    return (freq / samplerate).cumsum()

def pm_phase(samplerate, base_freq, mod_amount, mod):
    phase = numpy.full_like(mod, base_freq / samplerate).cumsum()
    return phase + mod_amount * mod

def dry_phase(samplerate, duration, base_freq):
    time = numpy.linspace(0, duration, int(samplerate * duration))
    return base_freq * time

def plot_error_to_base_freq():
    samplerate = 44100
    duration = 1
    base_freq = numpy.geomspace(20, 20000, 100)
    table_size = 1024
    lfo_freq = 2
    mod_amount = 70

    lfo = generate_wave(samplerate, duration, lfo_freq)

    cubic_error = []
    sinc_error = []
    for base in base_freq:
        print(f"\r{base}", end="")
        phase = fm_phase(samplerate, base, mod_amount, lfo)

        additive = additive_saw(samplerate, base, phase)
        cubic = yoshimi_saw_cubic(samplerate, base, phase, table_size)
        sinc = sinc_saw(samplerate, base, phase, table_size)

        cubic_error.append(mean_absolute_error(additive, cubic))
        sinc_error.append(mean_absolute_error(additive, sinc))

    pyplot.title("Sawtooth Wavetable Error")
    pyplot.plot(base_freq, cubic_error, label="cubic", color="red", alpha=0.5, lw=1)
    pyplot.plot(base_freq, sinc_error, label="sinc", color="blue", alpha=0.5, lw=1)
    pyplot.xscale("log")
    pyplot.xlabel("Base Frequency [Hz]")
    pyplot.yscale("log")
    pyplot.ylabel("Error (log scale)")
    pyplot.grid()
    pyplot.grid(which="minor", lw=1, alpha=0.1)
    pyplot.legend()
    pyplot.show()

# plot_error_to_base_freq()
# exit()

samplerate = 44100
duration = 1
base_freq = 1000
table_size = 1024
lfo_freq = 2
mod_amount = 70

lfo = generate_wave(samplerate, duration, lfo_freq)

# phase = dry_phase(samplerate, duration, base_freq)
phase = fm_phase(samplerate, base_freq, mod_amount, lfo)
# phase = pm_phase(samplerate, base_freq, mod_amount, lfo)

additive = additive_saw(samplerate, base_freq, phase)
naive = naive_saw(samplerate, phase, table_size)
linterp = yoshimi_saw_linterp(samplerate, base_freq, phase, table_size)
cubic = yoshimi_saw_cubic(samplerate, base_freq, phase, table_size)
sinc = sinc_saw(samplerate, base_freq, phase, table_size)

# errors = calc_error(additive, naive, linterp, cubic, sinc)
# pp = pprint.PrettyPrinter(indent=4)
# pp.pprint(errors)

# compare(additive, naive)
# compare(additive, linterp)
# compare(additive, cubic)
# compare(additive, sinc)

soundfile.write("snd/modAdditive.wav", normalize(additive), samplerate)
soundfile.write("snd/modNaive.wav", normalize(naive), samplerate)
soundfile.write("snd/modLinterp.wav", normalize(linterp), samplerate)
soundfile.write("snd/modCubic.wav", normalize(cubic), samplerate)
soundfile.write("snd/modSinc.wav", normalize(sinc), samplerate)
