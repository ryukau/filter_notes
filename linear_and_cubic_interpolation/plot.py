import numpy as np
import scipy.signal as signal
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt


def getIR(t, interpFunc, nPoint):
    """IR = Impulse Response"""
    sig = np.eye(nPoint)
    return [interpFunc(*sig[i], t) for i in range(nPoint)]


def plotResponse(interpFunc, nPoint):
    sampleRate = 2

    cmap = plt.get_cmap("plasma")
    fig, ax = plt.subplots(2, 2)
    for t in np.linspace(0, 1, 9):
        fir = getIR(t, interpFunc, nPoint)

        ω, response = signal.freqz(fir, fs=sampleRate, worN=2048)
        gain = 20 * np.log10(np.abs(response))
        phase = np.unwrap(np.angle(response))

        _, delay = signal.group_delay((fir, 1), fs=sampleRate, w=2048)

        ax[0][0].plot(ω, gain, color=cmap(t), label=f"{t:.3f}")
        ax[0][1].plot(ω, delay, color=cmap(t), label=f"{t:.3f}")
        ax[1][0].plot(ω, phase, color=cmap(t), label=f"{t:.3f}")
        ax[1][1].plot(fir, color=cmap(t), label=f"{t:.3f}")
    ax[0][0].set_ylabel("Gain [dB]")
    ax[0][0].set_xlabel(f"Normalized Frequency [rad/π]")
    ax[0][0].set_ylim((-12, 1))
    # ax[0][0].set_xscale("log")

    ax[0][1].set_ylabel("Group Delay [sample]")
    ax[0][1].set_xlabel(f"Normalized Frequency [rad/π]")
    ax[0][1].set_ylim((-0.5, 3.5))
    # ax[0][1].set_xscale("log")

    ax[1][0].set_ylabel("Phase [rad]")
    ax[1][0].set_xlabel(f"Normalized Frequency [rad/π]")
    # ax[1][0].set_xscale("log")

    ax[1][1].set_ylabel("Amplitude of FIR Coefficient")
    ax[1][1].set_xlabel(f"Time [sample]")
    ax[1][1].set_xticks(np.arange(nPoint))

    for row in ax:
        for axis in row:
            # axis.set_xlim((sampleRate / 8, sampleRate / 2))
            # axis.legend(ncol=3)
            axis.grid(which="both", color="#f0f0f0")
    fig.suptitle(f"{interpFunc.__name__}")
    fig.set_size_inches((10, 6))
    plt.tight_layout()
    # plt.show()
    plt.savefig(f"img/{interpFunc.__name__}.svg")
    plt.close()


def linear(y0, y1, t):
    return y0 + t * (y1 - y0)


def cubic_old(y0, y1, y2, y3, t):
    t2 = t * t
    c0 = y1 - y2
    c1 = (y2 - y0) * 0.5
    c2 = c0 + c1
    c3 = c0 + c2 + (y3 - y1) * 0.5
    return c3 * t * t2 - (c2 + c3) * t2 + c1 * t + y1


def cubic(y0, y1, y2, y3, t):
    """Using Horner's method."""
    c0 = y1 - y2
    c1 = (y2 - y0) * 0.5
    c2 = c0 + c1
    c3 = c0 + c2 + (y3 - y1) * 0.5
    return ((c3 * t - c2 - c3) * t + c1) * t + y1


def lagrange3(y0, y1, y2, y3, t):
    u = 1 + t
    d0 = y0 - y1
    d1 = d0 - (y1 - y2)
    d2 = d1 - ((y1 - y2) - (y2 - y3))
    return y0 - u * (d0 + (1 - u) / 2 * (d1 + (2 - u) / 3 * d2))


def pchip(y0, y1, y2, y3, t):
    m0 = y1 - y0
    m1 = y2 - y1
    m2 = y3 - y2

    dk0 = 0.0 if m0 * m1 <= 0 else 2 * (m0 * m1) / (m0 + m1)
    dk1 = 0.0 if m1 * m2 <= 0 else 2 * (m1 * m2) / (m1 + m2)

    t2 = t * t
    c0 = y1 - y2
    c1 = dk0
    c2 = c0 + c1
    c3 = c0 + c2 + dk1
    return c3 * t * t2 - (c2 + c3) * t2 + c1 * t + y1


def akima(y0, y1, y2, y3, y4, y5, t):
    """
    Interpolates y[3] and y[2].
    6 samples are required. y[5] is current sample, y[0] is 5 sample before.
    t is in [0, 1].

    - [akima/akima.py at master · cgohlke/akima · GitHub](https://github.com/cgohlke/akima/blob/master/akima/akima.py)
    - [scipy/_cubic.py at v1.8.0 · scipy/scipy · GitHub](https://github.com/scipy/scipy/blob/v1.8.0/scipy/interpolate/_cubic.py#L364-L461)
    """
    m0 = y1 - y0
    m1 = y2 - y1
    m2 = y3 - y2
    m3 = y4 - y3
    m4 = y5 - y4

    w2 = abs(m1 - m0)
    w3 = abs(m2 - m1)
    w4 = abs(m3 - m2)
    w5 = abs(m4 - m3)

    b2 = m1 if w2 + w4 == 0 else (w4 * m1 + w2 * m2) / (w2 + w4)
    b3 = m2 if w3 + w5 == 0 else (w5 * m2 + w3 * m3) / (w3 + w5)

    c2 = 3.0 * m2 - 2.0 * b2 - b3

    d2 = b2 + b3 - 2.0 * m2

    return ((t * d2 + c2) * t + b2) * t + y2


def uniformBSpline(y0, y1, y2, y3, t):
    t2 = t * t
    t3 = t2 * t

    b0 = 1 / 6 - t / 2 + t2 / 2 - t3 / 6
    b1 = 2 / 3 - t2 + t3 / 2
    b2 = 1 / 6 + t / 2 + t2 / 2 - t3 / 2
    b3 = t3 / 6
    return b0 * y0 + b1 * y1 + b2 * y2 + b3 * y3


def centripetalCatmullRom(y0, y1, y2, y3, t, alpha=0.1):
    if t <= 0:
        return y1
    if t >= 1:
        return y2

    #
    # abs is only applicable for 1D case. In 2D or higher, use dot product as following.
    # ```
    # d1 = y0 - y1
    # t1 = t0 + pow(dot(d1, d1), alpha / dimension)
    # # Same for t2 and t3.
    # ```
    #
    t0 = 0
    t1 = t0 + np.power(abs(y1 - y0), alpha)
    if t0 == t1:
        t1 += np.finfo(np.float64).eps
    t2 = t1 + np.power(abs(y2 - y1), alpha)
    if t1 == t2:
        t2 += t1 * np.finfo(np.float64).eps
    t3 = t2 + np.power(abs(y3 - y2), alpha)
    if t2 == t3:
        t3 += t2 * np.finfo(np.float64).eps

    # It must be t0 < t1 < t2 < t3 when reaching here. Otherwise, 0 division will occur.

    t = t1 + t * (t2 - t1)

    A1 = (t1 - t) / (t1 - t0) * y0 + (t - t0) / (t1 - t0) * y1
    A2 = (t2 - t) / (t2 - t1) * y1 + (t - t1) / (t2 - t1) * y2
    A3 = (t3 - t) / (t3 - t2) * y2 + (t - t2) / (t3 - t2) * y3
    B1 = (t2 - t) / (t2 - t0) * A1 + (t - t0) / (t2 - t0) * A2
    B2 = (t3 - t) / (t3 - t1) * A2 + (t - t1) / (t3 - t1) * A3
    return (t2 - t) / (t2 - t1) * B1 + (t - t1) / (t2 - t1) * B2


def centripetalCatmullRom01(y0, y1, y2, y3, t):
    return centripetalCatmullRom(y0, y1, y2, y3, t, 0.1)


def centripetalCatmullRom05(y0, y1, y2, y3, t):
    return centripetalCatmullRom(y0, y1, y2, y3, t, 0.5)


def centripetalCatmullRom09(y0, y1, y2, y3, t):
    return centripetalCatmullRom(y0, y1, y2, y3, t, 0.9)


def sinc4(y0, y1, y2, y3, t):
    s0 = np.sinc(t + 1) / np.e
    s1 = np.sinc(t)
    s2 = np.sinc(t - 1)
    s3 = np.sinc(t - 2) / np.e
    denom = s0 + s1 + s2 + s3
    return (y0 * s0 + y1 * s1 + y2 * s2 + y3 * s3) / denom


if __name__ == "__main__":
    plotResponse(linear, 2)
    plotResponse(cubic, 4)
    plotResponse(lagrange3, 4)
    plotResponse(pchip, 4)
    plotResponse(akima, 6)
    plotResponse(uniformBSpline, 4)
    # plotResponse(centripetalCatmullRom01, 4)
    plotResponse(centripetalCatmullRom05, 4)
    # plotResponse(centripetalCatmullRom09, 4)
    plotResponse(sinc4, 4)
    plt.show()
