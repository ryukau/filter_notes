import matplotlib.pyplot as pyplot
import numpy as np
import scipy.signal as signal


def frequencyShift(sig, normalizedFreq):
    analytic = signal.hilbert(sig)
    norm = np.abs(analytic)
    theta = np.angle(analytic)
    time = np.linspace(0, len(analytic), len(analytic))
    return norm * np.cos(theta + 2 * np.pi * normalizedFreq * time)


def toGroupDelay(response):
    gd = -np.diff(np.unwrap(np.angle(response))) * len(response) / np.pi
    # gd -= gd[0]
    return np.hstack(([gd[0]], gd)) - gd[0]


def getResponse(name, sig, fs):
    worN = np.geomspace(1, fs // 2, 2**10)
    freq, response = signal.freqz(sig, worN=worN, fs=fs)
    return {
        "name": name,
        "frequency": freq,
        "amp": 20 * np.log10(np.abs(response)),
        "phase": np.unwrap(np.angle(response)),
        "groupDelay": toGroupDelay(response),
    }


def plotImpulseResponse(samplerate, uprate, sosBandpass, sosLowpass):
    # sosBandpass = np.array(sosBandpass, dtype=np.float32)
    # sosLowpass = np.array(sosLowpass, dtype=np.float32)

    impulse = np.zeros(uprate)
    impulse[0] = 1

    ir = signal.sosfilt(sosBandpass, impulse)
    ir = frequencyShift(ir, -1 / 6)
    ir = 3 * signal.resample(ir, samplerate)  # Note the multiply by 3.

    lp = 6 * signal.sosfilt(sosLowpass, ir)
    mixed = lp + ir

    mixedResponse = getResponse("mixed", mixed, len(mixed))
    irResponse = getResponse("bp", ir, len(ir))
    lpResponse = getResponse("lp", lp, len(lp))
    freq = mixedResponse["frequency"]

    fig, ax = pyplot.subplots(3, 1)
    fig.set_size_inches((8, 9))

    for data in [mixedResponse, irResponse, lpResponse]:
        ax[0].set_title("Amp")
        ax[0].plot(freq, data["amp"], label=data["name"])
        ax[0].set_ylim((-20, 10))
        ax[1].set_title("Phase")
        ax[1].plot(freq, data["phase"], label=data["name"])
        ax[2].set_title("Group Delay")
        ax[2].plot(freq, data["groupDelay"], label=data["name"])
    for axis in ax:
        # axis.axvline(uprate * 1 / 6, lw=1, ls=":", color="red", alpha=0.5)
        # axis.axvline(uprate * 2 / 6, lw=1, ls=":", color="red", alpha=0.5)
        # axis.axvline(samplerate * 3 / 2, lw=1, ls=":", color="red", alpha=0.5)
        axis.axvline(samplerate / 2, lw=1, ls="-", color="red", alpha=0.5)
        axis.grid(which="both", color="#f0f0f0")
        axis.legend()
        axis.set_xscale("log")

    fig.tight_layout()
    pyplot.show()


if __name__ == "__main__":
    samplerate = 48000
    uprate = 3 * samplerate

    sosBandpass = signal.ellip(
        12,
        0.01,
        140,
        [1.02 * uprate / 6, 1.96 * uprate / 6],
        "bandpass",
        output="sos",
        fs=uprate,
    )

    sosLowpass = signal.butter(
        2,
        0.01 * samplerate / 2,
        "lowpass",
        output="sos",
        fs=samplerate,
    )

    plotImpulseResponse(samplerate, uprate, sosBandpass, sosLowpass)
