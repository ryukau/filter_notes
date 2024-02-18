import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def plotResponse(data: list, worN=8192, fs=48000):
    fig, ax = plt.subplots(4, 1)

    # worN = np.hstack([[0], np.geomspace(100, fs / 2, worN)])

    for fir in data:
        name = fir["name"]
        coefficient = fir["coefficient"]
        # print(name, coefficient)
        freq, resp = signal.freqz(coefficient, 1, worN=worN, fs=fs)
        freq, delay = signal.group_delay((coefficient, 1), w=worN, fs=fs)
        gain = 20 * np.log10(np.abs(resp))
        phase = np.unwrap(np.angle(resp))

        # print(f"delay: {delay[1]}")

        ax[0].plot(freq, gain, alpha=0.5, label=name)
        ax[1].plot(freq, gain, alpha=0.5, label=name)
        ax[2].plot(freq, phase, alpha=0.5, label=name)
        ax[3].plot(freq, delay, alpha=0.5, label=name)

    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim((-20, 3))
    ax[1].set_ylabel("Gain [dB]")
    ax[1].set_ylim((-200, 20))
    ax[2].set_ylabel("Phase [rad/sample]")
    ax[3].set_ylabel("Group Delay [sample]")
    ax[3].set_ylim([29, 33])
    ax[3].set_xlabel("Frequency [Hz]")
    fig.suptitle(f"Pre-filter (fs={fs}Hz)")

    for axis in ax:
        axis.axvline(fs / 2, color="black", ls="--", label="Nyquist")
        axis.axvline(
            fs * 3 / 8, color="black", alpha=0.2, ls="--", label=f"{fs * 3/8}Hz"
        )
        axis.grid(color="#f0f0f0", which="both")
        axis.legend(ncol=2, loc=3)
        # axis.set_xscale("log")
    fig.set_size_inches((8, 8))
    fig.tight_layout()
    plt.show()


def getInversedMaxUnderReadFIR(firLength=64):
    """`firLength` Should be even."""
    normalizedFreq = np.linspace(0, 0.5, firLength // 2 + 1, endpoint=True)
    invMaxUnderRead = np.cos(np.pi * normalizedFreq)
    fir = np.fft.irfft(invMaxUnderRead.astype(np.complex128))
    return np.roll(fir, len(fir) // 2 - 1)


def getBasicLimiterFIR():
    return signal.firwin(63, 18000, window=("dpss", 8), fs=48000)


def plotMaxUnderReadIdea():
    phase = np.linspace(0, 1, 1024)
    sin = np.sin(np.pi * phase)

    mid = 1 / 2
    truePeak = np.sin(np.pi * mid)

    fig = plt.figure()
    fig.set_size_inches((5, 4))
    plt.title(
        "Example of Maximum Under-read Condition\nAccording to ITU-R BS.1770 Definition (f=fs/4)"
    )
    plt.plot(2 * (phase - 0.25), sin, lw=1, color="black", label="Sine")
    plt.scatter(
        [mid],
        [truePeak],
        s=64,
        color="red",
        label="True-peak",
        zorder=100,
    )
    markerline, stemlines, baseline = plt.stem(
        [0, 1],
        [np.sin(np.pi * 0.25), np.sin(np.pi * 0.75)],
        basefmt=" ",
        label="Sample-peak",
    )
    plt.setp(markerline, markersize=8, color="black", zorder=100)
    plt.setp(stemlines, linestyle="--", alpha=0.5, color="black")
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.xticks(np.linspace(-0.5, 1.5, 5))
    plt.yticks([0, np.sin(np.pi * 0.75), 1])
    plt.grid(color="#e0e0e0")
    plt.legend(loc=8)
    plt.tight_layout()
    plt.show()


plotMaxUnderReadIdea()
plotResponse(
    [
        {"name": "MaxUnderRead", "coefficient": getInversedMaxUnderReadFIR()},
        {"name": "BasicLimiter", "coefficient": getBasicLimiterFIR()},
    ]
)
