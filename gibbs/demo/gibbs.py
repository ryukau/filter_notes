import numpy
import matplotlib.pyplot as pyplot
from scipy.signal import *
from scipy.special import eval_gegenbauer, gamma


def spectralLowpass(sig, numSeries):
    if numSeries >= len(sig):
        return sig
    spec = numpy.fft.rfft(sig)
    spec[numSeries + 1 : len(sig) - 1] = 0
    return numpy.fft.irfft(spec)


def additiveSaw(length, numSeries):
    return spectralLowpass(numpy.linspace(-1, 1, length), numSeries)


def additiveNoise(length, numSeries):
    return spectralLowpass(numpy.random.uniform(-1, 1, length), numSeries)


def analyticSignal(length, numSeries):
    sig = numpy.linspace(-1, 1, length)
    noise = numpy.zeros(length)
    order = 128
    rand = numpy.random.random(order)
    for i in range(order):
        noise += rand[i] * sig ** (i + 1)
    return spectralLowpass(noise, numSeries)


def applyFilter(source, filterFunc, power=1):
    filt = numpy.zeros(len(source))
    half = (len(source) - 1) / 2
    for i in range(len(filt)):
        eta = (i - half) / half
        filt[i] = filterFunc(eta) ** power
    filt /= numpy.sum(filt)

    return (convolve(source, filt, mode="full"), filt)


def plotGibbsSuppression(length, filterFunc, numSeries=16, power=1):
    source = additiveSaw(length, numSeries)
    suppressed, filt = applyFilter(source, filterFunc, power)
    source = numpy.append(source, numpy.zeros(length))
    plot(
        source,
        suppressed,
        filt,
        filterFunc.__name__ + " Filter",
        "result_" + filterFunc.__name__ + ".png",
    )


def fejer(eta):
    return 1 - eta


def lanczos(eta):
    return numpy.sinc(eta)


def raisedCosine(eta):
    return (1 + numpy.cos(numpy.pi * eta)) / 2


def sharpendRaisedCosine(eta):
    s = raisedCosine(eta)
    return s**4 * (35 - 84 * s + 70 * s**2 - 20 * s**3)


def exponential(eta, alpha=2**10, p=4):
    return numpy.exp(-alpha * eta**p)


def daubechies(eta):
    return 1 - 210 * eta**4 + 504 * eta**5 - 420 * eta**6 + 120 * eta**7


def plot(signal, result, params, param_name, file_name):
    fig = pyplot.figure(figsize=(6, 8))

    ax1 = fig.add_subplot(311)
    ax1.set_title("Signal")
    ax1.grid()
    ax1.plot(signal)

    ax2 = fig.add_subplot(312)
    ax2.set_title("Result")
    ax2.grid()
    ax2.plot(result)

    ax3 = fig.add_subplot(313)
    ax3.set_title(param_name)
    ax3.grid()
    ax3.plot(params)

    fig.tight_layout()
    fig.savefig(file_name)
    pyplot.close()


def plotGottleibShu(length, numSeries=16):
    source = additiveSaw(length, numSeries)
    suppressed, g_hat = gottliebShu(source, numSeries)
    plot(
        source,
        suppressed,
        g_hat,
        "Approximated Gegenbauer Coefficients",
        "result_gottleibShu.png",
    )


def plotGottleibShuNoise(length, signalFunc, numSeries=16):
    source = signalFunc(length, numSeries)
    suppressed, g_hat = gottliebShu(source, numSeries)
    regibbsed = spectralLowpass(suppressed, numSeries)

    plot(
        source,
        suppressed,
        regibbsed,
        "Filtered result",
        "result_gottleibShu_" + signalFunc.__name__ + ".png",
    )

    # Plot spectrum.
    source_spec = numpy.fft.rfft(source)[0 : numSeries + 1]
    regibbsed_spec = numpy.fft.rfft(regibbsed)[0 : numSeries + 1]
    source_power = numpy.abs(source_spec)
    regibbsed_power = numpy.abs(regibbsed_spec)
    source_phase = numpy.angle(source_spec)
    regibbsed_phase = numpy.angle(regibbsed_spec)

    fig = pyplot.figure(figsize=(8, 6))

    ax11 = fig.add_subplot(221)
    ax11.set_title("Signal power")
    ax11.grid()
    ax11.plot(source_power)

    ax21 = fig.add_subplot(222)
    ax21.set_title("Signal phase")
    ax21.grid()
    ax21.set_ylim(-numpy.pi, numpy.pi)
    ax21.plot(source_phase)

    ax12 = fig.add_subplot(223)
    ax12.set_title("Filtered result power")
    ax12.grid()
    ax12.plot(regibbsed_power)

    ax22 = fig.add_subplot(224)
    ax22.set_title("Filtered result phase")
    ax22.grid()
    ax22.set_ylim(-numpy.pi, numpy.pi)
    ax22.plot(regibbsed_phase)

    fig.tight_layout()
    fig.savefig("spectrum_comparison_" + signalFunc.__name__ + ".png")
    pyplot.close()


def gottliebShu(sig, N, gam=0.25):
    """
    :param sig: numpy.array input signal.
    :param N: Number of spectral partial sum or number of overtone.
    :param gam: Positive real constant.
    """
    lam = gam * N  # N が大きいと eval_gegenbauerでオーバーフローする。

    x = numpy.linspace(-1.0, 1.0, len(sig))
    k = numpy.arange(0, int(N / 2))
    ggnbr = eval_gegenbauer(numpy.tile(numpy.vstack(k), len(sig)), lam, x)
    ggnbr_t = numpy.transpose(ggnbr)

    h_k_lam = (
        numpy.sqrt(numpy.pi) * ggnbr_t[-1] * gamma(lam + 0.5) / gamma(lam) / (k + lam)
    )

    g_hat = (
        numpy.sum(numpy.power(1 - x * x, numpy.abs(lam - 0.5)) * sig * ggnbr, axis=1)
        / h_k_lam
    )

    return (numpy.sum(g_hat * ggnbr_t, axis=1), g_hat)


if __name__ == "__main__":
    plotGibbsSuppression(1024, fejer, 16, 1)
    plotGibbsSuppression(1024, lanczos, 16, 1)
    plotGibbsSuppression(1024, raisedCosine, 16, 1)
    plotGibbsSuppression(1024, sharpendRaisedCosine, 16, 1)
    plotGibbsSuppression(1024, exponential, 16, 1)
    plotGibbsSuppression(1024, daubechies, 16, 1)

    plotGottleibShu(1024, 16)
    plotGottleibShuNoise(1024, analyticSignal, 16)
    plotGottleibShuNoise(1024, additiveNoise, 16)
