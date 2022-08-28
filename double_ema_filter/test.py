import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

class EMA:
    def __init__(self):
        self.k = 1
        self.value = 0

    def process(self, x0):
        self.value += self.k * (x0 - self.value)
        return self.value

def testEma():
    sig = np.zeros(256)
    sig[0] = 1

    cutoff = 0.01
    ema = EMA()
    ema.k = cutoff

    emaOut = np.zeros_like(sig)
    for i, value in enumerate(sig):
        emaOut[i] = ema.process(value)

    def emaClosed(n, k):
        return k * (1 - k)**n

    closed = emaClosed(np.arange(len(sig)), cutoff)

    plt.plot(emaOut, lw=1, alpha=0.5, color="blue", label="emaOut")
    plt.plot(closed, lw=1, alpha=0.5, color="red", label="closed")
    plt.grid()
    plt.legend()
    plt.show()

class DoubleEMA:
    def __init__(self, k=1, value=0):
        self.k = k
        self.v0 = value
        self.v1 = value

    def reset(self, value=0):
        self.v0 = value
        self.v1 = value

    def process(self, x0):
        self.v0 += self.k * (x0 - self.v0)
        self.v1 += self.k * (self.v0 - self.v1)
        return self.v1

def testDoubleEma():
    sig = np.zeros(256)
    sig[0] = 1

    cutoff = 0.01
    ema = DoubleEMA()
    ema.k = cutoff

    emaOut = np.zeros_like(sig)
    for i, value in enumerate(sig):
        emaOut[i] = ema.process(value)

    def doubleEmaClosed(n, k):
        return k * k * (n + 1) * (1 - k)**n

    closed = doubleEmaClosed(np.arange(len(sig)), cutoff)

    # plt.plot(emaOut, lw=1, alpha=0.5, color="blue", label="emaOut")
    # plt.plot(closed, lw=1, alpha=0.5, color="red", label="closed")
    plt.plot(emaOut - closed, lw=1, alpha=0.5, color="red", label="diff")
    plt.grid()
    plt.legend()
    plt.show()

def testDoubleEmaStep():
    length = 256

    cutoff = 0.05
    ema = DoubleEMA()
    ema.k = cutoff

    emaOut = np.zeros(length)
    for i in range(length):
        emaOut[i] = ema.process(1)

    def doubleEmaClosedStep(n, k):
        value = 0
        for i in range(0, n + 1):
            value += k * k * (i + 1) * (1 - k)**i
        return value

    closed = np.array([doubleEmaClosedStep(n, cutoff) for n in np.arange(length)])

    plt.plot(emaOut, lw=1, alpha=0.5, color="blue", label="emaOut")
    plt.plot(closed, lw=1, alpha=0.5, color="red", label="closed")
    # plt.plot(emaOut - closed, lw=1, alpha=0.5, color="red", label="diff")
    plt.grid()
    plt.legend()
    plt.show()

def testDoubleEmaReverseStep():
    length = 256

    cutoff = 0.05
    ema = DoubleEMA()
    ema.k = cutoff
    ema.reset(1)

    emaOut = np.zeros(length)
    for i in range(length):
        emaOut[i] = ema.process(0)

    def doubleEmaClosedReverseStep(n, k, init, gain):
        value = 0
        for i in range(0, n + 1):
            value += k * k * (i + 1) * (1 - k)**i
        return init + gain * value

    closed = np.array(
        [doubleEmaClosedReverseStep(n, cutoff, 1, -1) for n in np.arange(length)])

    plt.plot(emaOut, lw=1, alpha=0.5, color="blue", label="emaOut")
    plt.plot(closed, lw=1, alpha=0.5, color="red", label="closed")
    # plt.plot(emaOut - closed, lw=1, alpha=0.5, color="red", label="diff")
    plt.grid()
    plt.legend()
    plt.show()

class DoubleEmaEnvelope():
    def __init__(self, k_A, k_D):
        self.attack = DoubleEMA(k_A, 0)
        self.decay = DoubleEMA(k_D, 1)

    def process(self):
        return self.attack.process(1) * self.decay.process(0)

def testDoubleEmaEnvelope():
    length = 256

    # A for attack, D for decay.
    k_A = 0.1
    k_D = 0.05

    emaA = DoubleEMA(k_A, 0)
    emaD = DoubleEMA(k_D, 1)
    emaOut = np.zeros(length)
    for i in range(length):
        emaOut[i] = emaA.process(1) * emaD.process(0)

    envelope = DoubleEmaEnvelope(k_A, k_D)
    envOut = np.zeros(length)
    for i in range(length):
        envOut[i] = envelope.process()

    def emaClosedSimplified(n, k, init, gain):
        value = (k * k * (n + 1) - k * n - 1) * (1 - k)**n + 1
        return init + gain * value

    closed = np.array([
        emaClosedSimplified(n, k_A, 0, 1) * emaClosedSimplified(n, k_D, 1, -1)
        for n in np.arange(length)
    ])

    plt.plot(envOut, lw=1, alpha=0.5, color="orange", label="envOut")
    plt.plot(emaOut, lw=1, alpha=0.5, color="blue", label="emaOut")
    plt.plot(closed, lw=1, alpha=0.5, color="red", label="closed")
    # plt.plot(emaOut - closed, lw=1, alpha=0.5, color="red", label="diff")
    # plt.ylim((-1e-15, 1e-15))
    plt.grid()
    plt.legend()
    plt.show()

def testClosed():
    def emaClosedSum(n, k, init, gain):
        value = 0
        for i in range(0, n + 1):
            value += k * k * (i + 1) * (1 - k)**i
        return init + gain * value

    def emaClosedSimplified(n, k, init, gain):
        value = (k * k * (n + 1) - k * n - 1) * (1 - k)**n + 1
        return init + gain * value

    length = 256
    cutoffA = 0.1
    cutoffD = 0.05

    closedSum = np.array([
        emaClosedSum(n, cutoffA, 0, 1) * emaClosedSum(n, cutoffD, 1, -1)
        for n in np.arange(length)
    ])
    closedSimplified = np.array([
        emaClosedSimplified(n, cutoffA, 0, 1) * emaClosedSimplified(n, cutoffD, 1, -1)
        for n in np.arange(length)
    ])

    plt.plot(closedSum, lw=1, alpha=0.5, color="blue", label="sum")
    plt.plot(closedSimplified, lw=1, alpha=0.5, color="red", label="simplified")
    # plt.plot(closedSum - closedSimplified, lw=1, alpha=0.5, color="red", label="diff")
    plt.grid()
    plt.legend()
    plt.show()

def testEmaDiff():
    def diff1(n, k, G):
        if k == 1:
            return 0
        return G * (k - 1) * (1 - k)**n * ((k * n + k + 1) * np.log(1 - k) + k)

    def diff2(n, k, G):
        if k == 1:
            return 0
        return -G * (1 - k)**(n + 1) * ((k * n + k + 1) * np.log(1 - k) + k)

    n = np.arange(256)
    k = 0.05
    gain = 1

    d1 = diff1(n, k, gain)
    d2 = diff2(n, k, gain)

    plt.plot(d1, lw=1, alpha=0.5, color="blue", label="d1")
    plt.plot(d2, lw=1, alpha=0.5, color="red", label="d2")
    # plt.plot(d1 - d2, lw=1, alpha=0.5, color="red", label="diff")
    plt.grid()
    plt.legend()
    plt.show()

def testDoubleEmaDiff():
    def diff1(n, k_A, k_D):
        return ((k_D**2 - k_D) * (-((k_A**2 * (n + 1) - k_A * n - 1) *
                                    (1 - k_A)**n + 1)) * (1 - k_D)**n -
                ((k_A**2 * (n + 1) - k_A * n - 1) * (1 - k_A)**n + 1) *
                (k_D**2 * (n + 1) - k_D * n - 1) * (1 - k_D)**n * np.log(1 - k_D) -
                (k_D**2 * (n + 1) - k_D * n - 1) * (1 - k_D)**n *
                ((k_A**2 - k_A) * (1 - k_A)**n + (k_A**2 * (n + 1) - k_A * n - 1) *
                 (1 - k_A)**n * np.log(1 - k_A)))

    def diff2(n, k_A, k_D):
        D_0 = (1 - k_D)**(n + 1)
        D_1 = k_D * n + k_D + 1
        D_2 = np.log(1 - k_D)
        A_0 = (1 - k_A)**(n + 1)
        A_1 = k_A * n + k_A + 1
        return D_0 * ((k_D + D_1 * D_2) * (1 - A_0 * A_1) - D_1 * A_0 *
                      (k_A + A_1 * np.log(1 - k_A)))

    def diff3(n, k_A, k_D):
        D_0 = (1 - k_D)**(n + 1)
        A_0 = (1 - k_A)**(n + 1)
        D_1 = k_D * n + k_D + 1
        A_1 = k_A * n + k_A + 1
        D_2 = np.log(1 - k_D)
        A_2 = np.log(1 - k_A)
        return D_0 * ((k_D + D_1 * D_2) * (1 - A_0 * A_1) - D_1 * A_0 * (k_A + A_1 * A_2))

    length = 256
    n = np.arange(length)
    k_A = 0.05
    k_D = 0.05

    envelope = DoubleEmaEnvelope(k_A, k_D)
    envOut = np.zeros(length)
    for i in range(length):
        envOut[i] = envelope.process()
    target = np.diff(envOut)

    d1 = diff1(n, k_A, k_D)
    d2 = diff2(n, k_A, k_D)
    d3 = diff3(n, k_A, k_D)

    # plt.plot(target, lw=1, alpha=0.5, color="black", label="target")
    plt.plot(d1, lw=1, alpha=0.5, color="blue", label="d1")
    plt.plot(d2, lw=1, alpha=0.5, color="red", label="d2")
    plt.plot(d3, lw=1, alpha=0.5, color="orange", label="d3")
    # plt.plot(d1 - d2, lw=1, alpha=0.5, color="red", label="diff")
    # plt.plot(d2 - d3, lw=1, alpha=0.5, color="red", label="diff")
    plt.grid()
    plt.legend()
    plt.show()

def testPeak():
    def doubleEmaEnvelopeD0(n, k_A, k_D):
        A = (1 - k_A)**(n + 1) * (k_A * n + k_A + 1)
        D = (1 - k_D)**(n + 1) * (k_D * n + k_D + 1)
        return (1 - A) * D

    def doubleEmaEnvelopeD1(n, k_A, k_D):
        D_0 = (1 - k_D)**(n + 1)
        A_0 = (1 - k_A)**(n + 1)
        D_1 = k_D * n + k_D + 1
        A_1 = k_A * n + k_A + 1
        D_2 = np.log(1 - k_D)
        A_2 = np.log(1 - k_A)
        return D_0 * ((k_D + D_1 * D_2) * (1 - A_0 * A_1) - D_1 * A_0 * (k_A + A_1 * A_2))

    length = 2048
    n = np.arange(length)
    k_A = 0.01
    k_D = 0.001

    envelope = DoubleEmaEnvelope(k_A, k_D)
    target = np.zeros(length)
    for i in range(length):
        target[i] = envelope.process()

    d0 = doubleEmaEnvelopeD0(n, k_A, k_D)
    d1 = doubleEmaEnvelopeD1(n, k_A, k_D)
    d2 = np.diff(d1)

    peak = np.argmax(d0)
    print(peak)

    # plt.plot(target, lw=10, alpha=0.5, color="yellow", label="target")
    # plt.plot(d0, lw=1, alpha=0.5, color="black", label="d0")
    plt.plot(d1, lw=1, alpha=0.5, color="blue", label="d1")
    # plt.plot(d2, lw=1, alpha=0.5, color="red", label="d2")
    plt.axvline(peak, lw=1, ls="--", color="black", alpha=0.5, label="peak")
    plt.grid(color="#f0f0f0")
    plt.legend()
    plt.show()

def samplesToKp(timeInSamples):
    y = 1 - np.cos(2 * np.pi / timeInSamples)
    return -y + np.sqrt(y * (y + 2))

def testFindPeak():
    def doubleEmaEnvelopeD0(n, k_A, k_D):
        A = (1 - k_A)**(n + 1) * (k_A * n + k_A + 1)
        D = (1 - k_D)**(n + 1) * (k_D * n + k_D + 1)
        return (1 - A) * D

    def doubleEmaEnvelopeD0Negative(n, k_A, k_D):
        A = (1 - k_A)**(n + 1) * (k_A * n + k_A + 1)
        D = (1 - k_D)**(n + 1) * (k_D * n + k_D + 1)
        return (A - 1) * D

    def doubleEmaEnvelopeD1(n, k_A, k_D):
        D_0 = (1 - k_D)**(n + 1)
        A_0 = (1 - k_A)**(n + 1)
        D_1 = k_D * n + k_D + 1
        A_1 = k_A * n + k_A + 1
        D_2 = np.log(1 - k_D)
        A_2 = np.log(1 - k_A)
        return D_0 * ((k_D + D_1 * D_2) * (1 - A_0 * A_1) - D_1 * A_0 * (k_A + A_1 * A_2))

    length = 48000
    n = np.arange(length)
    sampleA = 10 * 48000
    sampleD = 1 * 48000
    k_A = samplesToKp(sampleA)
    k_D = samplesToKp(sampleD)

    envelope = DoubleEmaEnvelope(k_A, k_D)
    target = np.zeros(length)
    for i in range(length):
        target[i] = envelope.process()

    d0 = doubleEmaEnvelopeD0(n, k_A, k_D)
    d1 = doubleEmaEnvelopeD1(n, k_A, k_D)

    #
    # Using bounds costs more iteration.
    #
    # Considering the computation cost of `doubleEmaEnvelopeD1`, `minimize_scalar`
    # with Brent method is faster than `toms748`.
    #
    result = optimize.minimize_scalar(
        doubleEmaEnvelopeD0Negative,
        args=(k_A, k_D),
        # bounds=(0, max([sampleA, sampleD])),
        # method="bounded",
    )
    peakTime = result.x
    peakValue = -result.fun

    # x0, result = optimize.toms748(
    #     doubleEmaEnvelopeD1,
    #     0,
    #     max([sampleA, sampleD]),
    #     args=(k_A, k_D),
    #     full_output=True,
    # )
    # peakTime = x0
    # peakValue = doubleEmaEnvelopeD0(peakTime, k_A, k_D)

    print(f"Iteration            : {result.nit}")
    # print(f"Iteration            : {result.iterations}")
    print(f"Value                : {peakValue}")
    print(f"Peak   Time [sample] : {peakTime}")
    print(f"Attack Time [sample] : {sampleA}")
    print(f"Decay  Time [sample] : {sampleD}")
    print(np.max(target), np.max(d0))

    target /= peakValue
    d0 /= peakValue

    plt.plot(target, lw=10, alpha=0.5, color="yellow", label="target")
    plt.plot(d0, lw=1, alpha=0.5, color="blue", label="d0")
    # plt.plot(d1, lw=1, alpha=0.5, color="red", label="d1")
    plt.axvline(peakTime, lw=1, ls="--", color="black", alpha=0.5, label="peakTime")
    plt.grid(color="#f0f0f0")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    # testEma()
    # testDoubleEma()
    # testDoubleEmaStep()
    # testDoubleEmaReverseStep()
    # testDoubleEmaEnvelope()
    # testClosed()
    # testEmaDiff()
    # testDoubleEmaDiff()
    # testPeak()
    testFindPeak()
