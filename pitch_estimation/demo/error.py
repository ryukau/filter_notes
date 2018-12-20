import itertools
import matplotlib.pyplot as pyplot
from matplotlib.colors import LogNorm
from pathlib import Path

from evaluation import *
from pitch import *

def hz_to_cent(frequency):
    return 1200 * numpy.log2(frequency / 440)

def mean_absolute_error(true_value, data_value):
    return numpy.nanmean(numpy.abs(true_value - data_value))

def plot_error_per_frequency(errors, frequencies, xname, xlabel, signal_type):
    cmap = pyplot.get_cmap("inferno")
    palette = [cmap(i) for i in numpy.linspace(0, 1, len(errors))]

    cmap_hatch = pyplot.get_cmap("binary")
    palette_hatch = [cmap_hatch(i) for i in numpy.linspace(0, 1, len(errors))]
    hatches = itertools.cycle(("-", "OOO", "/", "oo", None, "*", "\\\\", ".."))

    index = numpy.arange(len(frequencies))

    pyplot.figure(figsize=(19.20, 10.80), dpi=10)
    pyplot.title(f"Error for each {xname.title()} ({signal_type})")
    pyplot.xlabel(xlabel)
    pyplot.ylabel("Mean Absolute Error")
    pyplot.grid(axis="y", color="#eeeeee")
    # pyplot.yscale("log")
    for i, (method, error) in enumerate(errors.items()):
        pyplot.bar(
            index + i / (len(errors) + 1),
            height=error,
            width=1 / len(errors),
            hatch=next(hatches),
            color=palette[i],
            edgecolor=palette_hatch[i],
            zorder=100,
            label=method,
        )
    pyplot.xticks(
        index,
        [f"{freq:7.2f}" for freq in frequencies],
        rotation="vertical",
    )
    legend = pyplot.legend(loc=0)
    legend.set_zorder(101)
    pyplot.savefig(
        f"img/error_per_{xname.replace(' ', '_')}_{signal_type}.png", dpi=100)
    pyplot.close()

def format_method_name(name):
    name = name.replace("_", "-")
    name = name.upper()
    name = name.replace("-TYPE1", " Type I")
    name = name.replace("-TYPE2", " Type II")
    return name

def plot_error_sum(error_sum, signal_type):
    index = numpy.arange(len(error_sum))

    errors = [error for error in error_sum.values()]
    methods = [format_method_name(method) for method in error_sum.keys()]

    title_str = signal_type.replace("_", " ").title()

    pyplot.figure(figsize=(6.4, 4.8))
    pyplot.title(f"Mean Absolute Error to {title_str} ({signal_type})")
    pyplot.subplots_adjust(bottom=0.22, left=0.16)
    pyplot.ylabel("Mean Absolute Error")
    pyplot.grid(axis="y", color="#eeeeee")
    pyplot.bar(index, height=errors, color="red", zorder=100)
    pyplot.xticks(index, methods, rotation=30, ha="right")
    pyplot.savefig(f"img/error_{signal_type}.png")
    pyplot.close()

def plot_error_mesh(method, datas, signal_type, xtick, xlabel, ytick, ylabel):
    error_2d = extract_error_2d_from_method(datas, method)
    error_2d = numpy.log10(error_2d)
    error_mesh = numpy.ma.masked_where(~numpy.isfinite(error_2d), error_2d)

    pyplot.figure(figsize=(12.8, 7.2))
    pyplot.subplots_adjust(bottom=0.16, left=0.16)
    pyplot.title(
        f"Mean Absolute Error of {format_method_name(method)}  ({signal_type})")
    pyplot.xlabel(xlabel)
    pyplot.ylabel(ylabel)
    pyplot.xticks(
        numpy.arange(len(xtick)),
        [f"{xt:.3f}" for xt in xtick],
        rotation=90,
        ha="left",
        fontsize=8,
    )
    pyplot.yticks(
        numpy.arange(len(ytick)),
        [f"{yt:.3f}" for yt in ytick],
        va="bottom",
        fontsize=8,
    )
    pyplot.pcolormesh(error_mesh, cmap="plasma")
    colorbar = pyplot.colorbar()

    # colorbar のラベルを書き換える。対数スケールにはならない。
    # label = colorbar.ax.get_yaxis().get_ticklabels()
    # for lbl in label:
    #     lbl.set_text(f"$e^{{{lbl.get_text()}}}$")
    # colorbar.ax.get_yaxis().set_ticklabels(label)

    pyplot.savefig(f"img/error_{signal_type}_{method}.png")
    pyplot.close()

def extract_error_1d(datas):
    errors = {method: [] for method in datas[0].pitch.keys()}
    for data in datas:
        true_cent = hz_to_cent(data.frequency)
        for method, pitch in data.pitch.items():
            errors[method].append(
                mean_absolute_error(true_cent, hz_to_cent(pitch)))
    error_sum = {
        method: numpy.nanmean(error) for method, error in errors.items()
    }
    return (errors, error_sum)

def gather_error_2d(datas):
    errors = {method: [] for method in datas[0][0].pitch.keys()}
    for error_index, data in enumerate(datas):
        for error in errors.values():
            error.append(0)
        for i in range(0, len(data)):
            for method, pitch in data[i].pitch.items():
                mae = mean_absolute_error(
                    hz_to_cent(data[i].frequency), hz_to_cent(pitch))
                if numpy.isfinite(mae):
                    errors[method][error_index] += mae
    return errors

def extract_error_2d(datas):
    error_x = gather_error_2d(datas)
    error_y = gather_error_2d(datas.T)
    error_sum = {
        method: numpy.nanmean(error) for method, error in error_x.items()
    }
    return (error_x, error_y, error_sum)

def extract_error_2d_from_method(datas, method):
    error = numpy.zeros(datas.shape)
    for x in range(len(datas)):
        for y in range(len(datas[x])):
            error[x][y] = mean_absolute_error(
                hz_to_cent(datas[x][y].frequency),
                hz_to_cent(datas[x][y].pitch[method]),
            )
    return error

def plot_sin_wave(path):
    signal_type = path.stem
    datas = numpy.load(str(path))
    error_per_freq, error_sum = extract_error_1d(datas)
    frequencies = [data.frequency for data in datas]
    plot_error_per_frequency(
        error_per_freq,
        frequencies,
        "frequency",
        "Sine Wave Frequency[Hz]",
        signal_type,
    )
    plot_error_sum(error_sum, signal_type)

def plot_sin_with_noise(path):
    signal_type = path.stem

    datas = numpy.load(str(path))
    error_per_freq, error_per_ratio, error_sum = extract_error_2d(datas)

    frequencies = [data[0].frequency for data in datas]
    ratios = [data[0].misc["noise_ratio"] for data in datas.T]

    # mpm_nsd_type2 のエラーがとても大きいので別に分ける。
    # error_per_freq.pop("mpm_nsd_type2")
    # error_per_ratio.pop("mpm_nsd_type2")
    # error_sum.pop("mpm_nsd_type2")

    for func in pitch_functions:
        plot_error_mesh(
            func.__name__,
            datas,
            signal_type,
            ratios,
            "Noise Ratio",
            frequencies,
            "Frequency",
        )

    plot_error_per_frequency(
        error_per_freq,
        frequencies,
        "frequency",
        "Sine Wave Frequency[Hz]",
        signal_type,
    )
    plot_error_per_frequency(
        error_per_ratio,
        ratios,
        "noise ratio",
        "Noise Ratio",
        signal_type,
    )
    plot_error_sum(error_sum, signal_type)

def plot_sin_modulation(path):
    signal_type = path.stem

    datas = numpy.load(str(path))
    error_per_freq, error_per_ratio, error_sum = extract_error_2d(datas)

    frequencies = [data[0].frequency for data in datas]
    ratios = [data[0].misc["mod_freq"] for data in datas.T]

    for func in pitch_functions:
        plot_error_mesh(
            func.__name__,
            datas,
            signal_type,
            ratios,
            "Modulator Frequency [Hz]",
            frequencies,
            "Frequency",
        )

    plot_error_per_frequency(
        error_per_freq,
        frequencies,
        "carrier frequency",
        "Carrier Frequency [Hz]",
        signal_type,
    )
    plot_error_per_frequency(
        error_per_ratio,
        ratios,
        "modulator frequency",
        "Modulator Frequency [Hz]",
        signal_type,
    )
    plot_error_sum(error_sum, signal_type)

data_dir = Path("data")

# plot_sin_wave(data_dir / Path("sin_wave.npy"))
# plot_sin_with_noise(data_dir / Path("sin_with_noise.npy"))
plot_sin_modulation(data_dir / Path("sin_am.npy"))
plot_sin_modulation(data_dir / Path("sin_fm.npy"))
