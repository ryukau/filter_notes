import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

def samplesToKp(timeInSamples):
    y = 1 - np.cos(2 * np.pi / timeInSamples)
    return -y + np.sqrt(y * (y + 2))

def doubleEmaEnvelopeD0(n, k_A, k_D):
    A = (1 - k_A)**(n + 1) * (k_A * n + k_A + 1)
    D = (1 - k_D)**(n + 1) * (k_D * n + k_D + 1)
    return (1 - A) * D

def doubleEmaEnvelopeD0Negative(n, k_A, k_D):
    A = (1 - k_A)**(n + 1) * (k_A * n + k_A + 1)
    D = (1 - k_D)**(n + 1) * (k_D * n + k_D + 1)
    return (A - 1) * D

def plotNaiveEnvelope():
    length = 10000
    nData = 16

    n = np.arange(length)
    cmap = plt.get_cmap("plasma")

    k_D = samplesToKp(length)
    plt.figure(figsize=(8, 4))
    plt.title(f"AD Envelope Using Double EMA Filter with Varying k_A (k_D = {k_D:.3e})")
    for index, k_A in enumerate(samplesToKp(np.linspace(2, length, nData))):
        envelope = doubleEmaEnvelopeD0(n, k_A, k_D)
        plt.plot(envelope, color=cmap(index / nData), label=f"k_A={k_A:.3e}")
    plt.grid()
    plt.legend(ncol=2)
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.tight_layout()
    plt.ylim((-0.05, 1.05))

    k_A = samplesToKp(length)
    plt.figure(figsize=(8, 4))
    plt.title(f"AD Envelope Using Double EMA Filter with Varying k_D (k_A = {k_D:.3e})")
    for index, k_D in enumerate(samplesToKp(np.linspace(2, length, nData))):
        envelope = doubleEmaEnvelopeD0(n, k_A, k_D)
        plt.plot(envelope, color=cmap(index / nData), label=f"k_A={k_A:.3e}")
    plt.grid()
    plt.legend(ncol=2)
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.tight_layout()
    plt.ylim((-0.05, 1.05))

    plt.show()

def plotNormalizedEnvelope():
    length = 10000
    nData = 16

    n = np.arange(length)
    cmap = plt.get_cmap("plasma")

    k_D = samplesToKp(length)
    plt.figure(figsize=(8, 4))
    plt.title(
        f"Normalized AD Envelope Using Double EMA Filter with Varying k_A (k_D = {k_D:.3e})"
    )
    for index, k_A in enumerate(samplesToKp(np.linspace(2, length, nData))):
        result = optimize.minimize_scalar(doubleEmaEnvelopeD0Negative, args=(k_A, k_D))
        peakValue = -result.fun
        envelope = doubleEmaEnvelopeD0(n, k_A, k_D) / peakValue
        plt.plot(envelope, color=cmap(index / nData), label=f"k_A={k_A:.3e}")
    plt.grid()
    plt.legend(ncol=2)
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.tight_layout()
    plt.ylim((-0.05, 1.05))

    k_A = samplesToKp(length)
    plt.figure(figsize=(8, 4))
    plt.title(
        f"Normalized AD Envelope Using Double EMA Filter with Varying k_D (k_A = {k_D:.3e})"
    )
    for index, k_D in enumerate(samplesToKp(np.linspace(2, length, nData))):
        result = optimize.minimize_scalar(doubleEmaEnvelopeD0Negative, args=(k_A, k_D))
        peakValue = -result.fun
        envelope = doubleEmaEnvelopeD0(n, k_A, k_D) / peakValue
        plt.plot(envelope, color=cmap(index / nData), label=f"k_A={k_A:.3e}")
    plt.grid()
    plt.legend(ncol=2)
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.tight_layout()
    plt.ylim((-0.05, 1.05))

    plt.show()

def checkError():
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

    def getIterativeEnvelopePeak(length, k_A, k_D):
        ema_attack = DoubleEMA(k_A)
        ema_decay = DoubleEMA(k_D, 1)

        sig = np.empty(length)
        for i in range(length):
            sig[i] = ema_attack.process(1) * (ema_decay.process(0))
        index = np.argmax(sig)
        return (index, sig[index])

    def getClosedEnvelopePeak(n, k_A, k_D):
        result = optimize.minimize_scalar(doubleEmaEnvelopeD0Negative, args=(k_A, k_D))
        return (result.x, -result.fun)

    length = 100000

    diff = []
    n = np.arange(length)
    time_in_sample = np.linspace(10000, 100000, 11)
    for k_A in samplesToKp(time_in_sample):
        diff.append([])
        for k_D in samplesToKp(time_in_sample):
            iter_peak = getIterativeEnvelopePeak(length, k_A, k_D)
            closed_peak = getClosedEnvelopePeak(n, k_A, k_D)

            diff[-1].append(iter_peak[1] / closed_peak[1])

            if abs(iter_peak[0] - closed_peak[0]) >= 1:
                print(
                    f"k_A: {k_A}, k_D: {k_D}, x_iter: {iter_peak[0]}, x_closed: {closed_peak[0]}"
                )

    if np.any(np.array(diff) > 1):
        print("Error: Overshoot occured")

    cmap = plt.get_cmap("plasma")
    nData = len(time_in_sample)
    plt.title("Normalized Peak")
    for index in range(nData):
        plt.plot(
            time_in_sample,
            diff[index],
            color=cmap(index / nData),
            label=f"n_A={time_in_sample[index]:.0f}",
        )
    # plt.yscale("log")
    plt.xlabel("k_D")
    plt.ylabel("Amplitude")
    plt.grid()
    plt.legend()
    plt.show()

plotNaiveEnvelope()
plotNormalizedEnvelope()
# checkError()
