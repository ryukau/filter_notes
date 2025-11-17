import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import functools


def modifiedSinc(x, cutoff):
    if x == 0:
        return 2 * cutoff
    return np.sin(np.pi * 2 * cutoff * x) / (np.pi * x)


def lowpassFir(length: int, cutoff: float, fractionSample: float):
    mid = fractionSample - (length // 2 + length % 2)
    fir = np.zeros(length)
    for i in range(length):
        x = i + mid
        fir[i] = modifiedSinc(x, cutoff)
    return fir


def lowpassBiquadNaive(length: int, cutoff: float, fractionSample: float):
    mid = fractionSample - (length // 2 + length % 2)

    omega = 2 * np.pi * cutoff
    phi = mid * omega
    k = 2 * np.cos(omega)
    u1 = np.sin(phi - omega)
    u2 = np.sin(phi - 2 * omega)

    fir = np.zeros(length)
    for i in range(length):
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0

        x = i + mid
        if x == 0:
            fir[i] = 2 * cutoff
        else:
            fir[i] = u0 / (np.pi * x)
    return fir


def lowpassBiquadPade(length: int, cutoff: float, fractionSample: float):
    mid = fractionSample - (length // 2 + length % 2)

    omega = 2 * np.pi * cutoff
    phi = mid * omega
    k = 2 * np.cos(omega)
    u1 = np.sin(phi - omega)
    u2 = np.sin(phi - 2 * omega)

    fir = np.zeros(length)
    for i in range(length):
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0

        x = i + mid
        if abs(x) < 0.1:
            # Pade approximant of modified sinc around x = 0. Derived from following code:
            #
            # ```maxima
            # sinc: sin(2 * %pi * f_c * x) / (%pi * x);
            # pade(taylor(sinc, x, 0, 4), 2, 2);
            # ```
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = u0 / (np.pi * x)
    return fir


def lowpassBiquadTrigIdentity(length: int, cutoff: float, fractionSample: float):
    mid = fractionSample - (length // 2 + length % 2)

    omega = 2 * np.pi * cutoff
    phi = mid * omega

    cos_omega = np.cos(omega)
    sin_omega = np.sin(omega)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    k = 2 * cos_omega
    u1 = sin_phi * cos_omega - cos_phi * sin_omega

    cos_2omega = 2 * cos_omega * cos_omega - 1
    sin_2omega = 2 * sin_omega * cos_omega
    u2 = sin_phi * cos_2omega - cos_phi * sin_2omega

    fir = np.zeros(length)
    for i in range(length):
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0

        x = i + mid
        if abs(x) < 0.1:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = u0 / (np.pi * x)
    return fir


def lowpassReinsch(length: int, cutoff: float, fractionSample: float):
    mid = fractionSample - (length // 2 + length % 2)

    omega = 2 * np.pi * cutoff
    phi = mid * omega
    A = 2 * np.sin(omega / 2)
    u = np.sin(phi - omega)
    v = A * np.cos(phi - omega / 2)
    k = A * A

    fir = np.zeros(length)
    for i in range(length):
        u = u + v
        v = v - k * u

        x = i + mid
        if abs(x) < 0.01:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = u / (np.pi * x)
    return fir


def lowpassStableQuad(length: int, cutoff: float, fractionSample: float):
    mid = fractionSample - (length // 2 + length % 2)

    omega = 2 * np.pi * cutoff
    phi = mid * omega
    k1 = np.tan(omega / 2)
    k2 = np.sin(omega)
    u = np.cos(phi - omega)
    v = np.sin(phi - omega)

    fir = np.zeros(length)
    for i in range(length):
        w = u - k1 * v
        v += k2 * w
        u = w - k1 * v

        x = i + mid
        if abs(x) < 0.01:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = v / (np.pi * x)

    return fir


def lowpassLanczosSimple(
    length: int, cutoff: float, fractionSample: float, lanczos_a: float = 3
):
    mid = fractionSample - (length // 2 + length % 2)

    fir = np.zeros(length)
    for i in range(length):
        x = i + mid
        lanczos = modifiedSinc(x / lanczos_a, cutoff) / cutoff / 2
        fir[i] = modifiedSinc(x, cutoff) * lanczos
    return fir


def lowpassWelch(length: int, cutoff: float, fractionSample: float):
    mid = fractionSample - (length // 2 + length % 2)

    fir = np.zeros(length)
    M = max(1, (length - 1) / 2)
    for i in range(length):
        A = (i - M) / M
        window = 1 - A * A
        fir[i] = modifiedSinc(i + mid, cutoff) * window
    return fir


def lowpassLanczosBiquad(
    length: int, cutoff: float, fractionSample: float, lanczos_width: int = 8
):
    mid = fractionSample - (length // 2 + length % 2)

    o1_omega = 2 * np.pi * cutoff
    o1_phi = mid * o1_omega
    o1_k = 2 * np.cos(o1_omega)
    o1_u1 = np.sin(o1_phi - o1_omega)
    o1_u2 = np.sin(o1_phi - 2 * o1_omega)

    o2_omega = 2 * np.pi * cutoff / float(lanczos_width)
    o2_phi = mid * o2_omega
    o2_k = 2 * np.cos(o2_omega)
    o2_u1 = np.sin(o2_phi - o2_omega)
    o2_u2 = np.sin(o2_phi - 2 * o2_omega)

    fir = np.zeros(length)
    A = lanczos_width / (2 * cutoff * np.pi * np.pi)
    for i in range(length):
        o1_u0 = o1_k * o1_u1 - o1_u2
        o1_u2 = o1_u1
        o1_u1 = o1_u0

        o2_u0 = o2_k * o2_u1 - o2_u2
        o2_u2 = o2_u1
        o2_u1 = o2_u0

        x = i + mid
        if abs(x) < 0.1:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)

            fc_lanczos = cutoff / lanczos_width
            p = np.pi * fc_lanczos * x
            p *= p
            fir[i] *= (2 / 3) * fc_lanczos * (15 - 7 * p) / (5 + p)
        else:
            fir[i] = o1_u0 * o2_u0 * (A / (x * x))
    return fir


def lowpassBlackmanHarrisBiquad(length: int, cutoff: float, fractionSample: float):
    isEven = 1 - length % 2

    mid = fractionSample - (length // 2 + length % 2)

    o1_omega = 2 * np.pi * cutoff
    o1_phi = mid * o1_omega
    o1_k = 2 * np.cos(o1_omega)
    o1_u1 = np.sin(o1_phi - o1_omega)
    o1_u2 = np.sin(o1_phi - 2 * o1_omega)

    o2_omega = 2 * np.pi / float(length + isEven)
    o2_phi = np.pi / 2
    o2_k = 2 * np.cos(o2_omega)
    o2_u1 = np.sin(o2_phi - (1 - isEven) * o2_omega)
    o2_u2 = np.sin(o2_phi - (2 - isEven) * o2_omega)

    fir = np.zeros(length)
    for i in range(length):
        o1_u0 = o1_k * o1_u1 - o1_u2
        o1_u2 = o1_u1
        o1_u1 = o1_u0

        o2_u0 = o2_k * o2_u1 - o2_u2
        o2_u2 = o2_u1
        o2_u1 = o2_u0

        x = i + mid
        if abs(x) < 0.1:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = o1_u0 / (np.pi * x)

        window = 0.21747 + o2_u0 * (-0.45325 + o2_u0 * (0.28256 + o2_u0 * -0.04672))
        fir[i] *= window
    return fir


def lowpassBlackmanNuttallBiquad(length: int, cutoff: float, fractionSample: float):
    isEven = 1 - length % 2

    mid = fractionSample - (length // 2 + length % 2)

    o1_omega = 2 * np.pi * cutoff
    o1_phi = mid * o1_omega
    o1_k = 2 * np.cos(o1_omega)
    o1_u1 = np.sin(o1_phi - o1_omega)
    o1_u2 = np.sin(o1_phi - 2 * o1_omega)

    o2_omega = 2 * np.pi / float(length + isEven)
    o2_phi = np.pi / 2
    o2_k = 2 * np.cos(o2_omega)
    o2_u1 = np.sin(o2_phi - (1 - isEven) * o2_omega)
    o2_u2 = np.sin(o2_phi - (2 - isEven) * o2_omega)

    fir = np.zeros(length)
    for i in range(length):
        o1_u0 = o1_k * o1_u1 - o1_u2
        o1_u2 = o1_u1
        o1_u1 = o1_u0

        o2_u0 = o2_k * o2_u1 - o2_u2
        o2_u2 = o2_u1
        o2_u1 = o2_u0

        x = i + mid
        if abs(x) < 0.1:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = o1_u0 / (np.pi * x)

        window = 0.2269824 + o2_u0 * (
            -0.4572542 + o2_u0 * (0.273199 + o2_u0 * -0.0425644)
        )
        fir[i] *= window
    return fir


def lowpassHannBiquad(length: int, cutoff: float, fractionSample: float):
    isEven = 1 - length % 2

    mid = fractionSample - (length // 2 + length % 2)

    o1_omega = 2 * np.pi * cutoff
    o1_phi = mid * o1_omega
    o1_k = 2 * np.cos(o1_omega)
    o1_u1 = np.sin(o1_phi - o1_omega)
    o1_u2 = np.sin(o1_phi - 2 * o1_omega)

    o2_omega = np.pi / float(length)
    o2_k = 2 * np.cos(o2_omega)
    o2_u1 = np.sin((isEven - 1) * o2_omega)
    o2_u2 = np.sin((isEven - 2) * o2_omega)

    fir = np.zeros(length)
    for i in range(length):
        o1_u0 = o1_k * o1_u1 - o1_u2
        o1_u2 = o1_u1
        o1_u1 = o1_u0

        o2_u0 = o2_k * o2_u1 - o2_u2
        o2_u2 = o2_u1
        o2_u1 = o2_u0

        x = i + mid
        if abs(x) < 0.1:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = o1_u0 / (np.pi * x)

        window = o2_u0 * o2_u0
        fir[i] *= window
    return fir


def lowpassTriangleBiquad(length: int, cutoff: float, fractionSample: float):
    # Odd `length` uses bartlett window. Even `length` uses symmetric triangle.
    isEven = 1 - length % 2

    mid = fractionSample - (length // 2 + length % 2)

    omega = 2 * np.pi * cutoff
    phi = mid * omega
    k = 2 * np.cos(omega)
    u1 = np.sin(phi - omega)
    u2 = np.sin(phi - 2 * omega)

    fir = np.zeros(length)
    tri_N: int = max(2, length + isEven)
    inv_N: float = 1 / tri_N
    for i in range(length):
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0

        x = i + mid
        if abs(x) < 0.1:
            q = np.pi * cutoff * x
            q *= q
            fir[i] = (2 / 3) * cutoff * (15 - 7 * q) / (5 + q)
        else:
            fir[i] = u0 / (np.pi * x)

        window = 1 - abs(float(2 * (i + isEven) - tri_N) * inv_N)
        fir[i] *= window
    return fir


def lowpassReference(length: int, cutoff: float, fractionSample: float, windowType):
    mid = fractionSample - (length // 2 + length % 2)

    fir = np.zeros(length)
    for i in range(length):
        fir[i] = modifiedSinc(i + mid, cutoff)
    return fir * signal.get_window(windowType, length, fftbins=False)


def plotErrorAroundZero():
    length = 16
    cutoff = 0.1
    fraction = (1 + 0.3) * 2**-50

    accurate = lowpassFir(length, cutoff, fraction)
    naive = lowpassBiquadNaive(length, cutoff, fraction)
    trig = lowpassBiquadTrigIdentity(length, cutoff, fraction)
    pade = lowpassBiquadPade(length, cutoff, fraction)

    plt.plot(accurate, label="accurate", color="black", alpha=0.25, lw=4, ls="--")
    plt.plot(naive, label="naive", color="red", alpha=0.5, lw=1)
    plt.plot(trig, label="trig", color="blue", alpha=0.5, lw=1)
    plt.plot(pade, label="pade", color="orange", alpha=0.5, lw=2)
    plt.grid()
    plt.legend()
    plt.show()


def plotErrorVsFraction():
    length = 15
    cutoff = 0.5
    exponent = np.arange(-53, -1, 1, dtype=np.float64)
    fractionalDelay = np.linspace(0, 1, 16)

    cmap = plt.get_cmap("magma")
    for index, fd in enumerate(fractionalDelay):
        error = []
        for n in exponent:
            fraction = 2**n * (1 + fd)

            accurate = lowpassFir(length, cutoff, fraction)
            acc_max = np.max(np.abs(accurate))

            # fast = lowpassBiquadNaive(length, cutoff, fraction)
            # fast = lowpassBiquadTrigIdentity(length, cutoff, fraction)
            fast = lowpassBiquadPade(length, cutoff, fraction)
            fast_max = np.max(np.abs(fast))

            error.append(fast_max / acc_max)

        plt.plot(
            exponent,
            error,
            color=cmap(index / len(fractionalDelay)),
            alpha=0.5,
            lw=np.log1p(index) + 1,
            label=f"{fd:.3f}",
        )
    plt.xlabel("log2(fraction)")
    plt.ylabel("Relative Error (fast / accurate)")
    plt.legend(loc="center right")
    plt.grid()
    plt.show()


def compareInitialPhaseAccuracy():
    length = 256
    cutoff = 0.5
    mid = np.finfo(np.float64).eps - (length // 2 + length % 2)

    omega = 2 * np.pi * cutoff
    phi = mid * omega

    u1_naive = np.sin(phi - omega)
    u2_naive = np.sin(phi - 2 * omega)

    cos_omega = np.cos(omega)
    sin_omega = np.sin(omega)
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    u1_fixed = sin_phi * cos_omega - cos_phi * sin_omega
    cos_2omega = 2 * cos_omega * cos_omega - 1
    sin_2omega = 2 * sin_omega * cos_omega
    u2_fixed = sin_phi * cos_2omega - cos_phi * sin_2omega

    print(f"{u1_naive:+24}, {u2_naive:+24}")
    print(f"{u1_fixed:+24}, {u2_fixed:+24}")


def plotResponse(firList, cutoffList, nameList, worN=8192, fs=48000, title=""):
    fig, ax = plt.subplots(4, 1)
    if len(title) >= 1:
        fig.suptitle(title)

    # worN = np.hstack([[0], np.geomspace(0.1, fs / 2)])

    gdMedians = []
    cmap = plt.get_cmap("plasma")
    for idx, fir in enumerate(firList):
        freq, resp = signal.freqz(fir, 1, worN=worN, fs=fs)
        freq, delay = signal.group_delay((fir, 1), w=worN, fs=fs)
        gain = 20 * np.log10(np.abs(resp))
        phase = np.unwrap(np.angle(resp))

        cut = max(1, int(len(delay) * cutoffList[idx] / fs))
        # print(f"delay: {np.mean(delay[:cut])}")
        gdMedians.append(np.median(delay[:cut]))

        color = cmap(idx / len(firList))
        ax[0].axvline(cutoffList[idx], alpha=0.33, color="black", ls="--")
        ax[0].plot(freq, gain, alpha=0.66, lw=1, color=color)
        ax[1].plot(freq, phase, alpha=0.66, lw=1, color=color)
        ax[2].plot(freq, delay, alpha=0.66, lw=1, color=color)
        ax[3].plot(fir, alpha=0.66, lw=1, color=color, label=f"{nameList[idx]}")

    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim((-0.1, 0.1))
    # ax[0].set_ylim((-60, 6))
    ax[1].set_ylabel("Phase [rad/sample]")
    ax[2].set_ylabel("Delay [sample]")
    gdMid = np.median(gdMedians)
    gdRange = np.max(
        [
            1.5 * (gdMid - np.min(gdMedians)),
            1.5 * (np.max(gdMedians) - gdMid),
            1,
        ]
    )
    ax[2].set_ylim([gdMid - gdRange, gdMid + gdRange])
    ax[3].set_ylabel("FIR Amplitude")
    ax[3].legend(ncol=2)

    for axis in ax[:3]:
        # axis.set_xlim([100, 25000])
        axis.set_xscale("log")
        axis.axvline(fs / 2, color="black", ls="--")
        # axis.axvline(cutoffHz, color="black", ls="--")

    for axis in ax:
        axis.grid(color="#f0f0f0", which="both")
        # axis.legend(ncol=2)
        # axis.set_xscale("log")
    fig.set_size_inches((8, 8))
    fig.tight_layout()
    plt.show()


def compareResponseToReference():
    sampleRate = 1
    param = [256, 0.02, 0.6]  # [tap, cutoff, fraction]

    pairs = {
        # "boxcar": lowpassBiquadPade,
        # "hann": lowpassHannBiquad,
        # "blackmanharris": lowpassBlackmanHarrisBiquad,
        # "bartlett": lowpassTriangleBiquad,
        "welch": lowpassWelch,
    }

    name = ["ref.", "fast"]
    for key, value in pairs.items():
        fir = [lowpassReference(*param, key), value(*param)]
        plotResponse(
            fir,
            [param[1] * sampleRate] * len(fir),
            name,
            fs=sampleRate,
            title=f"{key}, cutoff={param[1]:.2f}, frac. delay={param[2]:.2f}",
        )


def compareResponseForDifferentWindows():
    sampleRate = 1
    param = [256, 0.025, 0.0001]  # [tap, cutoff, fraction]

    pairs = {
        # "boxcar": lowpassBiquadPade,
        # "hann": lowpassHannBiquad,
        "blackmanharris": lowpassBlackmanHarrisBiquad,
        # "blackmannuttall": lowpassBlackmanNuttallBiquad,
        # "bartlett": lowpassTriangleBiquad,
        # "welch": lowpassWelch,
        # "lanczos": functools.partial(lowpassLanczosBiquad, lanczos_width=4),
        # "hamming": functools.partial(lowpassReference, windowType="hamming"),
        # "hann": functools.partial(lowpassReference, windowType="hann"),
        # "blackman": functools.partial(lowpassReference, windowType="blackman"),
        # "nuttall": functools.partial(lowpassReference, windowType="nuttall"),
        "flattop": functools.partial(lowpassReference, windowType="flattop"),
    }

    name = list(pairs.keys())
    fir = [fn(*param) for fn in pairs.values()]
    plotResponse(
        fir,
        [param[1] * sampleRate] * len(pairs),
        name,
        fs=sampleRate,
        title=f"cutoff={param[1]}, frac. delay={param[2]}",
    )


if __name__ == "__main__":
    # plotErrorAroundZero()
    # plotErrorVsFraction()
    # compareInitialPhaseAccuracy()
    # compareResponseToReference()
    compareResponseForDifferentWindows()
