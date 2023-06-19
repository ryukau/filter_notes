import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import copy
from collections import deque


def toDecibel(response):
    return 20 * np.log10(np.maximum(np.abs(response), 1e-5))


def toPhase(response):
    return (np.angle(response)) / np.pi


def toGroupDelay(response):
    gd = -np.diff(np.unwrap(np.angle(response))) * len(response) / np.pi
    return gd - gd[0]


def biquadPair(cutoffNormalized, Q):
    """
    Returns `(lowpass_sos, highpass_sos)`. lowpass_sos and highpass_sos are biquad coefficients in sos (second order series) format.

    Reference:
    - https://www.w3.org/TR/audio-eq-cookbook/
    """

    ω0 = 2 * np.pi * cutoffNormalized
    cs = np.cos(ω0)
    α = np.sin(ω0) / (2 * Q)

    a0 = 1 + α
    a1 = (-2 * cs) / a0
    a2 = (1 - α) / a0

    lp_b0 = ((1 - cs) / 2) / a0
    lp_b1 = (1 - cs) / a0
    lp_b2 = ((1 - cs) / 2) / a0

    hp_b0 = ((1 + cs) / 2) / a0
    hp_b1 = -(1 + cs) / a0
    hp_b2 = ((1 + cs) / 2) / a0

    return (
        [lp_b0, lp_b1, lp_b2, 1, a1, a2],
        [hp_b0, hp_b1, hp_b2, 1, a1, a2],
    )


def biquadAllpass(cutoffNormalized, Q):
    ω0 = 2 * np.pi * cutoffNormalized
    cs = np.cos(ω0)
    α = np.sin(ω0) / (2 * Q)

    b2 = 1 + α
    b0 = (1 - α) / b2
    b1 = (-2 * cs) / b2

    return [b0, b1, 1, 1, b1, b0]


def onepoleAllpass(cutoffNormalized):
    t = np.tan(cutoffNormalized * np.pi)
    a = (t - 1) / (t + 1)
    return [a, 1, 0, 1, a, 0]


def testOnePoleAllpass():
    impulse = np.zeros(48000)
    impulse[0] = 1

    sampleRate = 2
    cutoff = 0.1
    ir = signal.sosfilt(onepoleAllpass(cutoff / sampleRate), impulse)
    freq, response = signal.freqz(ir, worN=2**16, fs=sampleRate)

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))
    ax[0].plot(freq, toDecibel(response), color="black")
    ax[1].plot(freq, toPhase(response), color="black")
    for axis in ax:
        axis.axvline(cutoff, color="black", ls="--", lw=1, label="$f_c$")
        axis.set_xscale("log")
        axis.grid(which="both", color="#f8f8f8")
        axis.legend()
    ax[0].set_title(f"1-pole Allpass Response ($f_s$=2, $f_c$={cutoff})")
    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim([-10, 10])
    ax[1].set_ylabel("Phase [rad/π]")
    ax[1].set_xlabel("Frequency [Hz]")
    fig.tight_layout()
    plt.show()


def butterAllpass(order, cutoffNormalized):
    sos = []
    if order % 2 == 1:
        sos.append(onepoleAllpass(cutoffNormalized))
        for k in range(1, order // 2 + 1):
            # Q is obtained from normalized Butterworth polynomials.
            # - https://en.wikipedia.org/wiki/Butterworth_filter
            Q = -2 * np.cos(np.pi * (2 * k + order - 1) / (2 * order))
            sos.append(biquadAllpass(cutoffNormalized, Q))
    else:
        K = order // 2
        for k in range(K):
            Q = 0.5 / np.sin((k + 0.5) * np.pi / K)
            sos.append(biquadAllpass(cutoffNormalized, Q))
    return sos


def testButterAllpass():
    sampleRate = 2
    order = 3
    cutoff = 0.5

    # Reference.
    lp_sos, hp_sos, _ = getLinkwitzRileySos(2 * order, cutoff / sampleRate)
    sig = np.zeros(48000)
    sig[0] = 1
    sign = 1 if order % 2 == 0 else -1
    target = signal.sosfilt(lp_sos, sig) + sign * signal.sosfilt(hp_sos, sig)
    freq, resp_ref = signal.freqz(target, worN=2**16, fs=sampleRate)

    # Test subject.
    sos = butterAllpass(order, cutoff / sampleRate)
    freq, resp_impl = signal.sosfreqz(sos, worN=2**16, fs=sampleRate)

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))
    ax[0].plot(freq, toDecibel(resp_ref), color="red", alpha=0.5, label="ref")
    ax[1].plot(freq, toPhase(resp_ref), color="red", alpha=0.5, label="ref")
    ax[0].plot(freq, toDecibel(resp_impl), color="blue", alpha=0.5, label="impl")
    ax[1].plot(freq, toPhase(resp_impl), color="blue", alpha=0.5, label="impl")
    for axis in ax:
        axis.axvline(cutoff, color="black", ls="--", lw=1, label="$f_c$")
        axis.set_xscale("log")
        axis.grid(which="both", color="#f8f8f8")
        axis.legend()
    ax[0].set_title(f"Allpass Test Response ($f_s$=2, $f_c$={cutoff})")
    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim([-10, 10])
    ax[1].set_ylabel("Phase [rad/π]")
    ax[1].set_xlabel("Frequency [Hz]")
    fig.tight_layout()
    plt.show()


class ComplexIIR:
    def __init__(self, pole, stage=4):
        pole = copy.copy(pole)

        self.stage = stage
        self.reset()
        self.pole1 = pole
        self.poles = []
        for i in range(stage - 1):
            pole *= pole
            self.poles.append(pole)

    def reset(self):
        self.x1 = 0
        self.delay = [
            deque([0 + 0j for _ in range(2**d)]) for d in range(1, self.stage)
        ]

    def process1PoleForward(self, x0):
        sig = x0 + self.pole1 * self.x1
        self.x1 = x0

        for k in range(len(self.poles)):
            self.delay[k].append(sig)
            sig = sig + self.poles[k] * self.delay[k].popleft()

        return sig

    def process1PoleReversed(self, x0):
        sig = self.pole1 * x0 + self.x1
        self.x1 = x0

        for k in range(len(self.poles)):
            self.delay[k].append(sig)
            sig = self.poles[k] * sig + self.delay[k].popleft()

        return sig

    def process2PoleForward(self, x0):
        sig = self.process1PoleForward(x0)
        return sig.real + (self.pole1.real / self.pole1.imag) * sig.imag

    def process2PoleReversed(self, x0):
        sig = self.process1PoleReversed(x0)
        return sig.real + (self.pole1.real / self.pole1.imag) * sig.imag


class LinearPhaseLinkwitzRiley4n:
    """
    Linear phase Linkwitz-Riley filter. Bilinear transform is used to obtain
    discrete transfer function.

    H(s) =         1 / ((1 - p * s) * (1 - conjugate(p) * s)).
    H(z) = (1 ± z)^2 / ((1 - p * z) * (1 - conjugate(p) * z)). <- this class.

    Filter becomes inaccurate for lower cutoff frequencies. More `stage` means
    better accuracy, in exchange of increased latency.

    Latency is `len(poles) * (2**stage - 1)` samples.
    """

    def __init__(self, gain, poles, stage, filterType="low"):
        poles = poles.copy()
        self.stage = stage
        self.forward = []
        self.reverse = []
        if filterType == "high":
            poles = np.conjugate(poles)
        else:
            assert filterType == "low", '`filterType` must be "low" or "high".'
        for pole in poles:
            self.reverse.append(ComplexIIR(pole, stage))
            self.forward.append(ComplexIIR(pole, stage))

        self.g1 = np.power(gain, 1 / len(poles))
        self.reset()

    def reset(self):
        N = len(self.reverse)

        self.u1 = np.zeros(N)
        self.u2 = np.zeros(N)
        self.v1 = np.zeros(N)
        self.v2 = np.zeros(N)

        for index in range(N):
            self.reverse[index].reset()
            self.forward[index].reset()

    def process(self, x0, a1):
        for idx in range(len(self.reverse)):
            u0 = self.reverse[idx].process2PoleReversed(x0 * self.g1)
            x0 = u0 + a1 * self.u1[idx] + self.u2[idx]
            self.u2[idx] = self.u1[idx]
            self.u1[idx] = u0

            v0 = self.forward[idx].process2PoleForward(x0 * self.g1)
            x0 = v0 + a1 * self.v1[idx] + self.v2[idx]
            self.v2[idx] = self.v1[idx]
            self.v1[idx] = v0

        return x0

    def processLowpass(self, x0):
        return self.process(x0, 2)

    def processHighpass(self, x0):
        return self.process(x0, -2)


class LinearPhaseLinkwitzRiley2:
    def __init__(self, gain, pole, stage, filterType="low"):
        if isinstance(pole, list) or isinstance(pole, np.ndarray):
            pole = pole[0]

        self.stage = stage
        if filterType == "high":
            pole = np.conjugate(pole)
        else:
            assert filterType == "low", '`filterType` must be "low" or "high".'
        self.forward = ComplexIIR(pole, stage)
        self.reverse = ComplexIIR(pole, stage)

        self.gain = gain
        self.reset()

    def reset(self):
        self.x1 = 0
        self.x2 = 0
        self.forward.reset()
        self.reverse.reset()

    def process(self, x0, a1):
        x0 = x0 * self.gain
        x0 = self.reverse.process1PoleReversed(x0)
        x0 = self.forward.process1PoleForward(x0)
        y0 = x0 + a1 * self.x1 + self.x2
        self.x2 = self.x1
        self.x1 = x0
        return y0

    def processLowpass(self, x0):
        return self.process(x0, 2)

    def processHighpass(self, x0):
        return self.process(x0, -2)


class LinearPhaseLinkwitzRileyApprox:
    """
    Approximate implementation of linear phase Linkwitz-Riley filter.

    This is basically a direct substitution of z to s in continuout time transfer function.

    H(s) = 1 / ((1 - p * s) * (1 - conjugate(p) * s)).
    H(z) = 1 / ((1 - p * z) * (1 - conjugate(p) * z)). <- this class.

    Highpass is subtraction of lowpass from input.
    """

    def __init__(self, gain, poles, stage, filterType="low"):
        self.stage = stage
        self.forward = []
        self.reverse = []
        if filterType == "high":
            poles = np.conjugate(poles)
        else:
            assert filterType == "low", '`filterType` must be "low" or "high".'
        for pole in poles:
            self.forward.append(ComplexIIR(pole, stage))
            self.reverse.append(ComplexIIR(pole, stage))

        self.g1 = np.power(gain, 1 / len(poles))
        self.g2 = np.power(gain * 4 ** len(poles), 1 / len(poles))
        self.reset()

    def reset(self):
        for index in range(len(self.reverse)):
            self.forward[index].reset()
            self.reverse[index].reset()

        self.hp_delay = deque(
            [0 + 0j for _ in range(len(self.reverse) * (2**self.stage - 1))]
        )

    def process(self, x0, a1):
        for index in range(len(self.reverse)):
            x0 = self.reverse[index].process2PoleReversed(x0 * self.g2)
            x0 = self.forward[index].process2PoleForward(x0 * self.g2)
        return x0

    def processLowpass(self, x0):
        return self.process(x0)

    def processHighpass(self, x0):
        self.hp_delay.append(x0)
        return self.hp_delay.popleft() - self.process(x0)


def getImpulseResponse(method: str, pole, stage):
    ir = np.zeros(2**16, dtype=np.complex128)
    ir[0] = 1
    forwardFilter = ComplexIIR(pole, stage)
    for i in range(len(ir)):
        ir[i] = getattr(forwardFilter, method)(ir[i])
    return ir


def testComplexIIROnePole():
    sampleRate = 2
    pole = 0.5 + 0.5j
    stage = 16

    assert pole.real != 0
    assert pole.imag != 0

    # Reference response.
    freq_ref, resp_ref = signal.freqz([1], [1, -pole], worN=2048, fs=sampleRate)

    # Impulse response of complex IIR.
    ir = getImpulseResponse("process1PoleForward", pole, stage)
    freq_fwd, resp_fwd = signal.freqz(ir, worN=2048, fs=sampleRate)

    # Impulse response of reversed IIR. `pole` is conjugated.
    ir = getImpulseResponse("process1PoleReversed", pole.conjugate(), stage)
    freq_rev, resp_rev = signal.freqz(ir, worN=2048, fs=sampleRate)

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))
    ref_label = "Reference"
    fwd_label = "Forward"
    rev_label = "Reversed"
    ax[0].plot(freq_ref, toDecibel(resp_ref), color="black", alpha=0.5, label=ref_label)
    ax[1].plot(freq_ref, toPhase(resp_ref), color="black", alpha=0.5, label=ref_label)
    ax[0].plot(freq_fwd, toDecibel(resp_fwd), color="red", alpha=0.5, label=fwd_label)
    ax[1].plot(freq_fwd, toPhase(resp_fwd), color="red", alpha=0.5, label=fwd_label)
    ax[0].plot(freq_rev, toDecibel(resp_rev), color="blue", alpha=0.5, label=rev_label)
    ax[1].plot(freq_rev, toPhase(resp_rev), color="blue", alpha=0.5, label=rev_label)
    for axis in ax:
        axis.set_xscale("log")
        axis.grid(which="both", color="#f8f8f8")
        axis.legend()
    ax[0].set_title(f"Complex 1-pole Responses, $f_s$={sampleRate}, $c$={pole}")
    ax[0].set_ylabel("Gain [dB]")
    ax[1].set_ylim([-1, 1])
    ax[1].set_ylabel("Phase [rad/π]")
    ax[1].set_xlabel("Frequency [Hz]")
    fig.tight_layout()
    plt.show()


def testComplexIIRTwoPole():
    sampleRate = 2
    pole = 0.15 + 0.57j
    stage = 16

    assert pole.real != 0
    assert pole.imag != 0

    # Reference response.
    freq_ref, resp_ref = signal.freqz(
        [1],
        [1, -(pole + pole.conjugate()), pole * pole.conjugate()],
        worN=2048,
        fs=sampleRate,
    )

    # Impulse response of complex IIR.
    ir = getImpulseResponse("process2PoleForward", pole, stage)
    freq_fwd, resp_fwd = signal.freqz(ir, worN=2048, fs=sampleRate)

    # Impulse response of reversed IIR. `pole` is conjugated.
    ir = getImpulseResponse("process2PoleReversed", pole, stage)
    freq_rev, resp_rev = signal.freqz(ir, worN=2048, fs=sampleRate)

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))
    ref_label = "Reference"
    fwd_label = "Forward"
    rev_label = "Reversed"
    ax[0].plot(freq_ref, toDecibel(resp_ref), color="black", alpha=0.5, label=ref_label)
    ax[1].plot(freq_ref, toPhase(resp_ref), color="black", alpha=0.5, label=ref_label)
    ax[0].plot(freq_fwd, toDecibel(resp_fwd), color="red", alpha=0.5, label=fwd_label)
    ax[1].plot(freq_fwd, toPhase(resp_fwd), color="red", alpha=0.5, label=fwd_label)
    ax[0].plot(freq_rev, toDecibel(resp_rev), color="blue", alpha=0.5, label=rev_label)
    ax[1].plot(freq_rev, toPhase(resp_rev), color="blue", alpha=0.5, label=rev_label)
    for axis in ax:
        axis.set_xscale("log")
        axis.grid(which="both", color="#f8f8f8")
        axis.legend()
    ax[0].set_title(f"Complex 2-pole Responses, $f_s$={sampleRate}, $c$={pole}")
    ax[0].set_ylabel("Gain [dB]")
    ax[1].set_ylim([-1, 1])
    ax[1].set_ylabel("Phase [rad/π]")
    ax[1].set_xlabel("Frequency [Hz]")
    fig.tight_layout()
    plt.show()


def plotLinkwitzRiley(
    order, sampleRate, cutoffHz, freq, response_lp, response_hp, response_sum
):
    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))

    ax[0].plot(freq, toDecibel(response_lp), color="red", alpha=0.5, label="LP")
    ax[0].plot(freq, toDecibel(response_hp), color="orange", alpha=0.5, label="HP")
    ax[0].plot(freq, toDecibel(response_sum), color="blue", alpha=0.5, label="Sum")

    ax[1].plot(freq, toPhase(response_lp), color="red", alpha=0.5, label="LP")
    ax[1].plot(freq, toPhase(response_hp), color="orange", alpha=0.5, label="HP")
    ax[1].plot(freq, toPhase(response_sum), color="blue", alpha=0.5, label="Sum")

    for axis in ax:
        axis.axvline(cutoffHz, color="black", ls="--", lw=1, label="$f_c$")
        axis.set_xscale("log")
        axis.grid(which="both", color="#f8f8f8")
        axis.legend()
    ax[0].set_title(f"Linkwitz-Riley {order} ($f_s$={sampleRate}, $f_c$={cutoffHz})")
    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim([-60, 10])
    ax[1].set_ylabel("Phase [rad/π]")
    ax[1].set_xlabel("Frequency [Hz]")
    fig.tight_layout()
    plt.show()


def plotLinkwitzRileySos(order, sampleRate, cutoffHz, lp_sos, hp_sos, worN=2048):
    freq, response_lp = signal.sosfreqz(lp_sos, worN=worN, fs=sampleRate)
    freq, response_hp = signal.sosfreqz(hp_sos, worN=worN, fs=sampleRate)
    response_sum = (
        response_lp + response_hp
        if (order // 2) % 2 == 0
        else response_lp - response_hp
    )
    plotLinkwitzRiley(
        order, sampleRate, cutoffHz, freq, response_lp, response_hp, response_sum
    )


def plotLinkwitzRileyIR(order, sampleRate, cutoffHz, lp_sos, hp_sos):
    impulse = np.zeros(int(sampleRate))
    impulse[0] = 1

    lp_out = signal.sosfilt(lp_sos, impulse)
    hp_out = signal.sosfilt(hp_sos, impulse)
    sum_out = lp_out + hp_out if (order // 2) % 2 == 0 else lp_out - hp_out

    freq = np.fft.rfftfreq(len(lp_out), 1 / sampleRate)
    response_lp = np.fft.rfft(lp_out)
    response_hp = np.fft.rfft(hp_out)
    response_sum = np.fft.rfft(sum_out)
    plotLinkwitzRiley(
        order, sampleRate, cutoffHz, freq, response_lp, response_hp, response_sum
    )


def getLinkwitzRileySos(order, normalizedCrossover):
    if order % 2 == 1:
        print("Odd order can't be specified for Linkwitz-Riley filter.")
        return

    lp_sos = []
    hp_sos = []
    ap_sos = []

    N = order // 2
    for n in range(N):
        Q = 0.5 / np.sin((n + 0.5) * np.pi / N)
        lp, hp = biquadPair(normalizedCrossover, Q)
        lp_sos.append(lp)
        hp_sos.append(hp)

    ap_sos = butterAllpass(N, normalizedCrossover)

    return (lp_sos, hp_sos, ap_sos)


def testLinkwitzRiley(order, sampleRate, cutoffHz):
    lp_sos, hp_sos, _ = getLinkwitzRileySos(order, cutoffHz / sampleRate)
    plotLinkwitzRileySos(order, sampleRate, cutoffHz, lp_sos, hp_sos)
    # plotLinkwitzRileyIR(order, sampleRate, cutoffHz, lp_sos, hp_sos)


def plotLinkwitzRileyCascade(bands, sum, order, sampleRate, cutoffsHz):
    freq = np.fft.rfftfreq(len(sum), 1 / sampleRate)
    plt.plot(freq, toDecibel(np.fft.rfft(sum)), color="black", label=f"sum")
    for idx, band in enumerate(bands):
        plt.plot(freq, toDecibel(np.fft.rfft(band)), alpha=0.5, label=f"band{idx}")

    for cut in cutoffsHz:
        plt.axvline(cut, color="black", ls="--")

    plt.title(f"LR{order}")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Gain [dB]")
    plt.xscale("log")
    plt.ylim([-60, 10])
    plt.grid(which="both", color="#f8f8f8")
    plt.legend()
    plt.tight_layout()
    plt.show()


def processLinkwitzRileyCascade(source, order, sampleRate, cutoffsHz):
    low = source.copy()

    bands = []
    for cut in [c / sampleRate for c in sorted(cutoffsHz, reverse=True)]:
        lp_sos, hp_sos, ap_sos = getLinkwitzRileySos(order, cut)
        bands = [signal.sosfilt(ap_sos, bd) for bd in bands]
        bands.append(signal.sosfilt(hp_sos, low))  # Higher band.
        low = signal.sosfilt(lp_sos, low)
    if order % 4 == 2:  # Invert sign for 2*odd order.
        bands = [-b for b in bands]
    bands.append(low)
    sum = np.sum(bands, axis=0)

    return (bands, sum)


def testLinkwitzRileyCascade(order, sampleRate, cutoffsHz):
    impulse = np.zeros(48000)
    impulse[0] = 1

    bands, sum = processLinkwitzRileyCascade(impulse, order, sampleRate, cutoffsHz)
    plotLinkwitzRileyCascade(bands, sum, order, sampleRate, cutoffsHz)


def convertPole(analogPole, cutoffHz, sampleRate=2):
    """
    Bilinear transform of all-pole analog prototype filter.

    Butterworth and Bessel are all-pole in analog. Use following routine to get
    prototype poles:

    - `scipy.signal.buttap`
    - `scipy.signal.besselap`

    To use the return values, utilize `zpk` (zero, pole gain) representation of
    scipy.signal.

    z = Empty because filter is all-pole.
    p = `poles`.
    k = `lp_gain` or `hp_gain`.

    `poles` represent 2-pole filter where both poles are complex conjugate to
    each other.

    H_cc(z) = 1 / ((1 - c * z^-1) * (1 - conjugate(c) * z^-1))
    """

    cutoffRadian = 2 * np.pi * cutoffHz / sampleRate
    poles = np.zeros_like(analogPole, dtype=np.complex128)
    co = np.zeros((len(poles), 2))
    lp_gain = 1
    hp_gain = 1
    for idx in range(len(poles)):
        pole = cutoffRadian * analogPole[idx]
        pole = (2 + pole) / (2 - pole)  # Bilinear transform.
        poles[idx] = pole
        co[idx][0] = -2 * pole.real  # Denominator coefficient for z^-1.
        co[idx][1] = pole.real * pole.real + pole.imag * pole.imag  # Same for z^-2.
        lp_gain *= (1 + co[idx][0] + co[idx][1]) / 4
        hp_gain *= (1 - co[idx][0] + co[idx][1]) / 4

    lp_sos = []
    hp_sos = []
    for coef in co:
        lp_sos.append([1, 2, 1, 1, coef[0], coef[1]])
        hp_sos.append([1, -2, 1, 1, coef[0], coef[1]])
    for i in range(3):
        lp_sos[0][i] *= lp_gain
        hp_sos[0][i] *= hp_gain

    return (poles, lp_gain, hp_gain, lp_sos, hp_sos)


def testConvertPole():
    order = 4
    sampleRate = 48000
    cutoffHz = 1000
    worN = 2048

    # Reference response.
    lp_sos_ref = signal.butter(order, cutoffHz, "low", output="sos", fs=sampleRate)
    hp_sos_ref = signal.butter(order, cutoffHz, "high", output="sos", fs=sampleRate)
    freq, resp_lp = signal.sosfreqz(lp_sos_ref, worN=worN, fs=sampleRate)
    freq, resp_hp = signal.sosfreqz(hp_sos_ref, worN=worN, fs=sampleRate)

    # Test response.
    _, p, _ = signal.buttap(order)
    p = list(sorted(p, key=lambda x: np.angle(x)))  # Just in case.
    p = p[: len(p) // 2]
    poles, lp_gain, hp_gain, lp_sos_blt, hp_sos_blt = convertPole(
        p, cutoffHz, sampleRate
    )
    freq, resp_blt_lp = signal.sosfreqz(lp_sos_blt, worN=worN, fs=sampleRate)
    freq, resp_blt_hp = signal.sosfreqz(hp_sos_blt, worN=worN, fs=sampleRate)

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))

    ax[0].plot(freq, toDecibel(resp_lp), color="red", alpha=0.5, label="Ref. LP")
    ax[0].plot(freq, toDecibel(resp_hp), color="orange", alpha=0.5, label="Ref. HP")
    ax[0].plot(freq, toDecibel(resp_blt_lp), color="yellow", alpha=0.5, label="BLT LP")
    ax[0].plot(freq, toDecibel(resp_blt_hp), color="purple", alpha=0.5, label="BLT HP")

    ax[1].plot(freq, toPhase(resp_lp), color="red", alpha=0.5, label="Ref. LP")
    ax[1].plot(freq, toPhase(resp_hp), color="orange", alpha=0.5, label="Ref. HP")
    ax[1].plot(freq, toPhase(resp_blt_lp), color="yellow", alpha=0.5, label="BLT LP")
    ax[1].plot(freq, toPhase(resp_blt_hp), color="purple", alpha=0.5, label="BLT HP")

    for axis in ax:
        axis.axvline(cutoffHz, color="black", ls="--", lw=1, label="$f_c$")
        axis.set_xscale("log")
        axis.grid(which="both", color="#f8f8f8")
        axis.legend()
    ax[0].set_title(f"Butterworth {order} ($f_s$={sampleRate}, $f_c$={cutoffHz})")
    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim([-60, 10])
    ax[1].set_ylabel("Phase [rad/π]")
    ax[1].set_xlabel("Frequency [Hz]")
    fig.tight_layout()
    plt.show()


def processLinearPhaseLinkwitzRiley(order, sampleRate, cutoffHz, stage, filterType):
    flt = []

    _, p, _ = signal.buttap(order // 2)
    p = list(sorted(p, key=lambda x: np.angle(x)))  # Just in case.

    isLowpass = filterType == "low"
    isOdd = len(p) % 2 == 1
    if isOdd:
        poles, lp_gain, hp_gain, _, _ = convertPole([p[0]], cutoffHz, sampleRate)
        gain = lp_gain if isLowpass else hp_gain
        flt.append(LinearPhaseLinkwitzRiley2(gain, poles, stage, filterType))
    if len(p) >= 2:
        if isOdd:
            p = p[1 : len(p) // 2 + 1]
        else:
            p = p[: len(p) // 2]
        poles, lp_gain, hp_gain, _, _ = convertPole(p, cutoffHz, sampleRate)
        gain = lp_gain if isLowpass else hp_gain
        flt.append(LinearPhaseLinkwitzRiley4n(gain, poles, stage, filterType))

    method = "processLowpass" if isLowpass else "processHighpass"
    ir = np.zeros(2**16, np.complex128)
    ir[0] = 1
    for i in range(len(ir)):
        for f in flt:
            ir[i] = getattr(f, method)(ir[i])
    return ir


def testLinearPhaseLinkwitzRiley(order, sampleRate, cutoffHz):
    worN = 2**16
    stage = 10

    sum_sign = 1 if (order // 2) % 2 == 0 else -1

    # Reference response.
    lp_sos, hp_sos, _ = getLinkwitzRileySos(order, cutoffHz / sampleRate)
    freq, resp_lp = signal.sosfreqz(lp_sos, worN=worN, fs=sampleRate)
    freq, resp_hp = signal.sosfreqz(hp_sos, worN=worN, fs=sampleRate)
    resp_sum = resp_lp + resp_hp * sum_sign

    # Test response.
    ir_lin_lp = processLinearPhaseLinkwitzRiley(
        order, sampleRate, cutoffHz, stage, "low"
    )
    ir_lin_hp = processLinearPhaseLinkwitzRiley(
        order, sampleRate, cutoffHz, stage, "high"
    )

    ir_lin_sum = ir_lin_lp + ir_lin_hp * sum_sign
    # plt.plot(ir_lin_hp)
    # plt.show()
    # exit()

    freq, resp_lin_lp = signal.freqz(ir_lin_lp, worN=worN, fs=sampleRate)
    freq, resp_lin_hp = signal.freqz(ir_lin_hp, worN=worN, fs=sampleRate)
    freq, resp_lin_sum = signal.freqz(ir_lin_sum, worN=worN, fs=sampleRate)

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 6))

    ax[0].plot(freq, toDecibel(resp_lp), color="red", alpha=0.5, label="Ref. LP")
    ax[0].plot(freq, toDecibel(resp_hp), color="red", alpha=0.5, label="Ref. HP")
    ax[0].plot(freq, toDecibel(resp_sum), color="red", alpha=0.5, label="Ref. Sum")
    ax[0].plot(freq, toDecibel(resp_lin_lp), color="blue", alpha=0.5, label="Lin. LP")
    ax[0].plot(freq, toDecibel(resp_lin_hp), color="blue", alpha=0.5, label="Lin. HP")
    ax[0].plot(freq, toDecibel(resp_lin_sum), color="blue", alpha=0.5, label="Lin. Sum")

    def plotGd(response, color, label):
        ax[1].plot(
            freq[1:], toGroupDelay(response), color=color, alpha=0.5, label=label
        )

    # plotGd(resp_lp, "red", "Ref. LP")
    # plotGd(resp_hp, "red", "Ref. HP")
    plotGd(resp_sum, "red", "Ref. Sum")
    # plotGd(resp_lin_lp, "blue", "Lin. LP")
    # plotGd(resp_lin_hp, "blue", "Lin. HP")
    plotGd(resp_lin_sum, "blue", "Lin. Sum")

    for axis in ax:
        axis.axvline(cutoffHz, color="black", ls="--", lw=1, label="$f_c$")
        axis.set_xscale("log")
        axis.grid(which="both", color="#f8f8f8")
        axis.legend()
    ax[0].set_title(
        f"Linear Phase Linkwitz-Riley {order} ($f_s$={sampleRate}, $f_c$={cutoffHz})"
    )
    ax[0].set_ylabel("Gain [dB]")
    ax[0].set_ylim([-60, 10])
    ax[1].set_ylabel("Δ Group Delay [sample]")
    ax[1].set_xlabel("Frequency [Hz]")
    ax[1].set_ylim([-10, 10])
    fig.tight_layout()
    plt.show()


def testPolarity(order, sampleRate, cutoffsHz):
    freq = np.geomspace(np.pi / 2**13, np.pi, 2**16)
    source = np.sin(freq.cumsum())

    _, sum = processLinkwitzRileyCascade(source, order, sampleRate, cutoffsHz)

    plt.plot(source)
    plt.plot(sum)
    plt.show()


def processSvf(sig, normalizedFreq, Q):
    """
    TPT biquad filter described in "The Art of VA Filter Design".
    """
    out_lp = np.empty_like(sig)
    out_hp = np.empty_like(sig)

    s1 = 0
    s2 = 0
    g = np.tan(normalizedFreq * np.pi)
    k = 1 / Q
    for idx, v0 in enumerate(sig):
        v1 = (s1 + g * (v0 - s2)) / (1 + g * (g + k))
        v2 = s2 + g * v1
        s1 = 2 * v1 - s1
        s2 = 2 * v2 - s2

        out_lp[idx] = v2
        out_hp[idx] = v0 - k * v1 - v2

    return (out_lp, out_hp)


if __name__ == "__main__":
    # testOnePoleAllpass()
    # testButterAllpass()
    # testComplexIIROnePole()
    # testComplexIIRTwoPole()
    # testConvertPole()

    sampleRate = 48000
    cutoffHz = 1000

    # for order in range(2, 16 + 1, 2):
    #     testLinkwitzRiley(order,sampleRate, cutoffHz)
    # testLinkwitzRiley(2, sampleRate, cutoffHz)

    # for order in range(2, 16 + 1, 2):
    #     testLinkwitzRileyCascade(order, sampleRate, [100, 500, 1000])
    # testLinkwitzRileyCascade(6, sampleRate, [100, 500, 1000])

    for order in range(2, 16 + 1, 2):
        testLinearPhaseLinkwitzRiley(order, sampleRate, cutoffHz)
    # testLinearPhaseLinkwitzRiley(8, sampleRate, cutoffHz)

    # testPolarity(4, sampleRate, [100, 500, 1000])
