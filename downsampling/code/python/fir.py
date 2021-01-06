import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

def ampToDecibel(z):
    return 20 * np.log10(np.abs(z))

def decimationFir(length, oversample, cutoff, nyquist=0.5, fs=48000, plot=True):
    bands = np.hstack((
        [0, cutoff],
        [nyquist, oversample / 2],
    ))
    desired = (1, 0)
    fir = signal.remez(length, bands, desired, fs=oversample)

    if plot:
        freq, h0 = signal.freqz(fir, worN=2048, fs=oversample * fs)
        plt.plot(freq, ampToDecibel(h0), label=f"{oversample}x")
        # plt.xscale("log")
        plt.grid(which="both")
        plt.legend()
        plt.show()

    return fir

def plotSos(sos, oversample, fs, worN=2048):
    fig, ax = plt.subplots(3, 1)
    fig.set_size_inches(6, 8)
    fig.set_tight_layout(True)

    ax[0].set_title("Amplitude Response")
    ax[0].set_ylabel("Amplitude [dB]")
    ax[1].set_title("Phase Response")
    ax[1].set_ylabel("Phase [rad]")
    ax[2].set_title("Group Delay")
    ax[2].set_ylabel("Group Delay [sample]")
    ax[2].set_xlabel("Frequency [rad/sample]")

    cmap = plt.get_cmap("plasma")

    freq, h0 = signal.sosfreqz(sos, worN=worN, fs=oversample * fs)

    b, a = signal.sos2tf(sos)
    _, gd = signal.group_delay((b, a), w=worN)

    ax[0].plot(freq, ampToDecibel(h0), alpha=0.75, color="black")
    ax[0].set_ylim((-140, 5))
    ax[1].plot(freq, np.unwrap(np.angle(h0)), alpha=0.75, color="black")
    ax[2].plot(freq, gd, alpha=0.75, color="black")

    for axis in ax:
        # axis.set_xscale("log")
        axis.grid(which="both")

    plt.show()

def decimationElliptic(
    order,
    oversample,
    cutoff,
    ripple=0.01,
    attenuation=120,
    fs=48000,
    plot=True,
):
    sos = signal.ellip(
        order,
        ripple,
        attenuation,
        cutoff,
        "low",
        output="sos",
        fs=oversample,
    )

    if plot:
        plotSos(sos, oversample, fs)

    return sos

def decimationButter(
    order,
    oversample,
    cutoff,
    fs=48000,
    plot=True,
):
    sos = signal.butter(
        order,
        cutoff,
        "low",
        output="sos",
        fs=oversample,
    )

    if plot:
        plotSos(sos, oversample, fs)

    return sos

def formatFir(fir, oversample):
    text = ""
    for phase in range(oversample):
        part = fir[phase::oversample]

        text += "{"
        for value in part:
            text += f"Sample({value}),"
        text = text[:-1]
        text += "},\n"
    print(text)

def formatSos(sos):
    text = ""
    for section in sos:
        text += "{"
        for idx, value in enumerate(section):
            if idx == 3:
                continue
            text += f"Sample({value}),"
        text = text[:-1]
        text += "},\n"
    print(text)

oversample = 4

fir = decimationFir(108, oversample, 0.375, 0.5, fs=48000)
formatFir(fir, oversample)

# sos = decimationElliptic(12, oversample, 0.4, ripple=0.01, attenuation=100, fs=1)
# sos = decimationButter(16, oversample, 0.4, fs=1)
# formatSos(sos)
