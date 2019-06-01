import matplotlib.pyplot as pyplot
import multiprocessing
from evaluation import *
from pitch import *

def parabolic_interpolation_plot(array, x):
    x_result = None
    if x < 1:
        x_result = x if array[x] <= array[x + 1] else x + 1
    elif x >= len(array) - 1:
        x_result = x if array[x] <= array[x - 1] else x - 1
    else:
        denom = array[x + 1] + array[x - 1] - 2 * array[x]
        delta = array[x - 1] - array[x + 1]
        if denom == 0:
            return (x, array[x])
        return (
            x + delta / (2 * denom),
            array[x] - delta * delta / (8 * denom),
        )
    return (x_result, array[x_result])

### MPM ###
def get_key_maxima(diff):
    start_index = 0
    while start_index < len(diff) and diff[start_index] > 0:
        start_index += 1
    isNegative = True
    max_index = 0
    key_maxima = []
    for i in range(start_index, len(diff)):
        if isNegative:
            if diff[i] < 0:
                continue
            max_index = i
            isNegative = False
        if diff[i] < 0:
            isNegative = True
            key_maxima.append((max_index, diff[max_index]))
        if diff[i] > diff[max_index]:
            max_index = i
    return key_maxima  # [(index, value), ...]

def pick_peak(diff):
    key_maxima = get_key_maxima(diff)
    if len(key_maxima) <= 0:
        return None
    maximum = max(key_maxima, key=lambda x: x[1])
    threshold = MPM_K * maximum[1]
    for index, value in key_maxima:
        if value >= threshold:
            return index
    warnings.warn("MPM_K may be greater than 1.")
    return maximum[0]

### Plot ###
def generate_sin(duration, samplerate, frequency):
    length = int(duration * samplerate)
    phase = numpy.linspace(0, 2 * numpy.pi * frequency * duration, length)
    return numpy.sin(phase)

def plot_autocorr():
    duration = 0.1
    frequency = 66
    samplerate = 16000

    sig = generate_sin(duration, samplerate, frequency)
    acf1 = autocorrelation_type1(sig)
    acf2 = autocorrelation_type2(sig)

    time = [t / samplerate for t in range(len(sig))]

    fig = pyplot.figure(figsize=(6.4, 6.4))

    ax1 = fig.add_subplot(211)
    ax1.set_title("Signal")
    ax1.grid(color="#eeeeee")
    ax1.plot(time, sig, lw=1, color="black")

    ax2 = fig.add_subplot(212)
    ax2.set_title("Autocorrelation")
    ax2.grid(color="#eeeeee")
    ax2.plot(time, acf1, lw=1, alpha=0.7, color="red", label="Type I")
    ax2.plot(time, acf2, lw=1, alpha=0.7, color="blue", label="Type II")
    ax2.legend()

    fig.tight_layout()
    pyplot.savefig("img/autocorr.png")
    pyplot.close()

def plot_diff():
    duration = 0.1
    frequency = 66
    samplerate = 16000

    sig = generate_sin(duration, samplerate, frequency)
    diff1 = difference_type1(sig)
    diff2 = difference_type2(sig)
    cmnd1 = cumulative_mean_normalized_difference(numpy.copy(diff1))
    cmnd2 = cumulative_mean_normalized_difference(numpy.copy(diff2))

    time = [t / samplerate for t in range(len(sig))]

    fig = pyplot.figure(figsize=(6.4, 6.4))

    ax = fig.add_subplot(211)
    ax.set_title("Difference of Autocorrelation")
    ax.grid(color="#eeeeee")
    ax.plot(time, diff1, lw=1, alpha=0.7, color="red", label="Type I")
    ax.plot(time, diff2, lw=1, alpha=0.7, color="blue", label="Type II")
    ax.legend()

    ax = fig.add_subplot(212)
    ax.set_title("Cumulative Mean Normalized Difference")
    ax.grid(color="#eeeeee")
    ax.plot(time, cmnd1, lw=1, alpha=0.7, color="red", label="Type I")
    ax.plot(time, cmnd2, lw=1, alpha=0.7, color="blue", label="Type II")
    ax.legend()

    fig.tight_layout()
    pyplot.savefig("img/diff.png")
    pyplot.close()

def plot_yin():
    duration = 0.1
    frequency = 66
    samplerate = 16000

    sig = generate_sin(duration, samplerate, frequency)
    cmnd1 = cumulative_mean_normalized_difference(difference_type1(sig))
    cmnd2 = cumulative_mean_normalized_difference(difference_type2(sig))

    local_min_cmnd1 = absolute_threshold(cmnd1)
    local_min_cmnd2 = absolute_threshold(cmnd2)

    est_i1, est_min1 = parabolic_interpolation_plot(cmnd1, local_min_cmnd1)
    est_i2, est_min2 = parabolic_interpolation_plot(cmnd2, local_min_cmnd2)

    print(f"type I: {samplerate / est_i1}Hz")
    print(f"type II: {samplerate / est_i2}Hz")

    time = [t / samplerate for t in range(len(sig))]

    fig = pyplot.figure(figsize=(6.4, 3.6))

    ax = fig.add_subplot(111)
    ax.set_title("YIN Pitch Estimation")
    ax.grid(color="#eeeeee")
    ax.plot(time, cmnd1, lw=1, alpha=0.7, color="red", label="Type I CMND")
    ax.plot(time, cmnd2, lw=1, alpha=0.7, color="blue", label="Type II CMND")
    ax.axhline(
        y=YIN_THRESHOLD,
        lw=1,
        ls="-.",
        color="black",
        alpha=0.5,
        label=f"Threshold",
    )
    ax.axvline(
        x=1 / frequency,
        lw=1,
        ls="--",
        color="black",
        alpha=0.5,
        label=f"True Period",
    )
    ax.scatter(
        est_i1 / samplerate,
        est_min1,
        s=48,
        color="red",
        zorder=100,
        marker=6,
        label=f"Type I Local Minimum")
    ax.scatter(
        est_i2 / samplerate,
        est_min2,
        s=48,
        color="blue",
        zorder=100,
        marker=6,
        label=f"Type II Local Minimum")
    ax.legend(loc=1)

    fig.tight_layout()
    pyplot.savefig("img/yin.png")
    pyplot.close()

def plot_mpm():
    duration = 0.1
    frequency = 66
    samplerate = 16000

    sig = generate_sin(duration, samplerate, frequency)
    nsd1 = normalized_square_difference_type1(numpy.copy(sig))
    nsd2 = normalized_square_difference_type2(numpy.copy(sig))

    index1 = estimate_period(nsd1)
    index2 = estimate_period(nsd2)

    est_i1, est_max1 = parabolic_interpolation_plot(nsd1, index1)
    est_i2, est_max2 = parabolic_interpolation_plot(nsd2, index2)

    print(f"type I: {samplerate / est_i1}Hz")
    print(f"type II: {samplerate / est_i2}Hz")

    time = [t / samplerate for t in range(len(sig))]

    fig = pyplot.figure(figsize=(6.4, 6.4))

    ax1 = fig.add_subplot(211)
    ax1.set_title("Signal")
    ax1.grid(color="#eeeeee")
    ax1.plot(time, sig, lw=1, color="black")

    ax = fig.add_subplot(212)
    ax.set_title("MPM Pitch Estimation")
    ax.grid(color="#eeeeee")
    ax.plot(time, nsd1, lw=1, alpha=0.7, color="red", label="Type I NSD")
    ax.plot(time, nsd2, lw=1, alpha=0.7, color="blue", label="Type II NSD")
    ax.axhline(
        y=MPM_K * numpy.max(nsd1),
        lw=1,
        ls="-.",
        color="red",
        alpha=0.2,
        label=f"Type I Threshold",
    )
    ax.axhline(
        y=MPM_K * numpy.max(nsd2),
        lw=1,
        ls="-.",
        color="blue",
        alpha=0.2,
        label=f"Type II Threshold",
    )
    ax.axvline(
        x=1 / frequency,
        lw=1,
        ls="--",
        color="black",
        alpha=0.5,
        label=f"True Period")
    ax.scatter(
        est_i1 / samplerate,
        est_max1,
        s=48,
        color="red",
        zorder=100,
        marker=7,
        label=f"Type I Local Maximum")
    ax.scatter(
        est_i2 / samplerate,
        est_max2,
        s=48,
        color="blue",
        zorder=100,
        marker=7,
        label=f"Type II Local Maximum")
    ax.legend(loc=1)

    fig.tight_layout()
    pyplot.savefig("img/mpm.png")
    pyplot.close()

def plot_yin_nsd():
    duration = 0.1
    frequency = 66
    samplerate = 16000

    sig = generate_sin(duration, samplerate, frequency)
    nsd1 = invert_nsd(normalized_square_difference_type1(sig))
    nsd2 = invert_nsd(normalized_square_difference_type2(sig))

    threshold = 0
    local_min_nsd1 = absolute_threshold(nsd1, threshold)
    local_min_nsd2 = absolute_threshold(nsd2, threshold)

    est_i1, est_min1 = parabolic_interpolation_plot(nsd1, local_min_nsd1)
    est_i2, est_min2 = parabolic_interpolation_plot(nsd2, local_min_nsd2)

    print(f"type I: {samplerate / est_i1}Hz")
    print(f"type II: {samplerate / est_i2}Hz")

    time = [t / samplerate for t in range(len(sig))]

    fig = pyplot.figure(figsize=(6.4, 3.6))

    ax = fig.add_subplot(111)
    ax.set_title("YIN-NSD Pitch Estimation")
    ax.grid(color="#eeeeee")
    ax.plot(time, nsd1, lw=1, alpha=0.7, color="red", label="Type I NSD")
    ax.plot(time, nsd2, lw=1, alpha=0.7, color="blue", label="Type II NSD")
    ax.axhline(
        y=threshold,
        lw=1,
        ls="-.",
        color="black",
        alpha=0.5,
        label=f"Threshold",
    )
    ax.axvline(
        x=1 / frequency,
        lw=1,
        ls="--",
        color="black",
        alpha=0.5,
        label=f"True Period",
    )
    ax.scatter(
        est_i1 / samplerate,
        est_min1,
        s=48,
        color="red",
        zorder=100,
        marker=6,
        label=f"Type I Local Minimum")
    ax.scatter(
        est_i2 / samplerate,
        est_min2,
        s=48,
        color="blue",
        zorder=100,
        marker=6,
        label=f"Type II Local Minimum")
    ax.legend(loc=4)

    fig.tight_layout()
    pyplot.savefig("img/yin_nsd.png")
    pyplot.close()

def plot_mpm_cmnd():
    duration = 0.1
    frequency = 66
    samplerate = 16000

    sig = generate_sin(duration, samplerate, frequency)
    cmnd1 = cumulative_mean_normalized_difference(difference_type1(sig))
    cmnd1 = numpy.max(cmnd1) / 2 - cmnd1
    cmnd2 = cumulative_mean_normalized_difference(difference_type2(sig))
    cmnd2 = numpy.max(cmnd2) / 2 - cmnd2

    local_min_cmnd1 = estimate_period(cmnd1)
    local_min_cmnd2 = estimate_period(cmnd2)

    est_i1, est_min1 = parabolic_interpolation_plot(cmnd1, local_min_cmnd1)
    est_i2, est_min2 = parabolic_interpolation_plot(cmnd2, local_min_cmnd2)

    print(f"type I: {samplerate / est_i1}Hz")
    print(f"type II: {samplerate / est_i2}Hz")

    time = [t / samplerate for t in range(len(sig))]

    fig = pyplot.figure(figsize=(6.4, 3.6))

    ax = fig.add_subplot(111)
    ax.set_title("MPM-CMND Pitch Estimation")
    ax.grid(color="#eeeeee")
    ax.plot(time, cmnd1, lw=1, alpha=0.7, color="red", label="Type I CMND")
    ax.plot(time, cmnd2, lw=1, alpha=0.7, color="blue", label="Type II CMND")
    ax.axhline(
        y=MPM_K * numpy.max(cmnd1),
        lw=1,
        ls="-.",
        color="red",
        alpha=0.2,
        label=f"Type I Threshold",
    )
    ax.axhline(
        y=MPM_K * numpy.max(cmnd2),
        lw=1,
        ls="-.",
        color="blue",
        alpha=0.2,
        label=f"Type II Threshold",
    )
    ax.axvline(
        x=1 / frequency,
        lw=1,
        ls="--",
        color="black",
        alpha=0.5,
        label=f"True Period",
    )
    ax.scatter(
        est_i1 / samplerate,
        est_min1,
        s=48,
        color="red",
        zorder=100,
        marker=7,
        label=f"Type I Local Maximum")
    ax.scatter(
        est_i2 / samplerate,
        est_min2,
        s=48,
        color="blue",
        zorder=100,
        marker=7,
        label=f"Type II Local Maximum")
    ax.legend(loc=4)

    fig.tight_layout()
    pyplot.savefig("img/mpm_cmnd.png")
    pyplot.close()

def plot_pitch(data, title):
    duration = len(data.signal) / data.samplerate
    time = [i / data.samplerate for i in range(len(data.signal))]

    cmap = pyplot.get_cmap("nipy_spectral")

    fig = pyplot.figure(figsize=(12.8, 7.2))

    ax = fig.add_subplot(211)
    ax.set_title(title)
    ax.set_ylabel("Amplitude")
    ax.set_xlim([0, duration])
    ax.grid(color="#eeeeee")
    ax.plot(time, data.signal, lw=1, color="black")

    ax = fig.add_subplot(212)
    ax.set_title("Estiamted Pitch")
    ax.set_xlabel("Time [s]", ha="left", x=0)
    ax.set_ylabel("Pitch [Hz]")
    ax.set_xlim([0, duration])
    ax.set_yscale("log")
    ax.set_ylim([10, 10000])
    ax.grid(color="#eeeeee", which="both")
    palette = [
        cmap(i) for i in numpy.linspace(1, 0, len(data.pitch), endpoint=False)
    ]
    len_pitch = len(list(data.pitch.values())[0])
    len_method = len(list(data.pitch.keys()))
    time = [i * data.winstep for i in range(len_pitch)]
    for i, (name, pitch) in enumerate(data.pitch.items()):
        ax.plot(
            time,
            pitch,
            lw=(len_method - i) * 1.5,
            alpha=0.5,
            color=palette[i],
            label=name)
    ax.axhline(
        y=data.frequency, lw=1, ls="--", color="black", label=f"Mod Frequency")
    ax.axhline(
        y=data.misc["mod_freq"],
        lw=1,
        ls="-.",
        color="black",
        label=f"Car Frequency")
    # ax.legend(loc=1)
    ax.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=5)

    fig.tight_layout()
    pyplot.savefig(f"img/pitch/{title}.png")
    # pyplot.show()
    pyplot.close()

def job_sin_mod(data, car_freq, mod_freq, index):
    title = f"{index:08d} AM Sin - Carrier {car_freq:.3f}Hz, Modulator {mod_freq:.3f}Hz"
    plot_pitch(data, title)

def plot_some_sin_mod():
    datas = numpy.load("data/sin_am.npy")
    index = 0
    for data in datas:
        indices = index + numpy.arange(len(data))
        index += len(data)
        args = [(d, d.frequency, d.misc["mod_freq"], i)
                for d, i in zip(data, indices)]
        with multiprocessing.Pool() as pool:
            pool.starmap(job_sin_mod, args)

def plot_sin_wave():
    datas = numpy.load("data/sin_wave.npy")
    for data in datas:
        plot_pitch(data)

# plot_autocorr()
plot_diff()
# plot_yin()
# plot_mpm()
# plot_yin_nsd()
# plot_mpm_cmnd()
# plot_sin_wave()
# plot_some_sin_mod()
