import numpy as np
from scipy.signal import lfilter, freqz
import matplotlib.pyplot as plt


def ema_alpha(cutoff_normalized):
    sn = np.sin(np.pi * cutoff_normalized)
    result = 2 * sn / (np.sqrt(sn * sn + 1) + sn)
    return result


class Lp3:
    def __init__(self):
        self.reset()

    def reset(self):
        self.s0 = 0.0
        self.s1 = 0.0
        self.s2 = 0.0
        self.x1 = 0.0

    def process(self, x0, A, B):
        self.s0 = A * self.s1 + B * self.s0
        self.s1 -= self.s0 + x0 - self.x1
        self.s2 -= self.s1 * A / (1.0 - B)
        self.x1 = x0
        return self.s2


class ChamberlinLp3:
    def __init__(self):
        self.reset()

    def reset(self):
        self.lp = 0
        self.bp = 0

    def processA(self, x0, A, B):
        f = np.sqrt(A)
        q = (1 - B) / f

        hp = x0 - self.lp - q * self.bp
        self.bp += f * hp
        self.lp += f * self.bp
        return self.lp

    def process(self, x0, A, B):
        f = np.sqrt(A)
        q = (1 - B) / f

        hp = x0 - self.lp - q * self.bp
        self.bp += f * hp
        self.lp += f * self.bp

        return self.lp + B * self.bp / q


class ChamberlinSVF:
    def __init__(self, bounded=True):
        self.lp = 0
        self.bp = 0
        self.cutoff_max = 1 / 6 if bounded else 0.4997

    def process(self, x0, cutoffNormalized, Q):
        """
        cutoffNormalized in [0, 0.5).
        """
        cut = np.clip(cutoffNormalized, 0, self.cutoff_max)
        f = 2 * np.sin(np.pi * cut)

        Q = max(Q, np.finfo(np.float64).eps)  # 0 除算を避ける。
        Q = max(Q, 2 * f / ((2 - f) * (2 + f)))

        hp = x0 - self.lp - self.bp / Q
        self.bp += f * hp
        self.lp += f * self.bp
        return self.lp


def plot_equivalence():
    alpha = 0.01
    beta = 0.9
    k = alpha / (1.0 - beta)

    np.random.seed(45287)
    signal_len = 1000
    input_sig = np.random.randn(signal_len)
    input_sig[0] = 10.0

    lp3 = Lp3()
    output_lp3 = np.array([lp3.process(x, alpha, beta) for x in input_sig])

    chamberlin = ChamberlinLp3()
    output_chamberlin = np.array(
        [chamberlin.process(x, alpha, beta) for x in input_sig]
    )

    b = [k, -k * beta, 0]
    a = [1, -(1 + beta - alpha), beta]
    output_tf = lfilter(b, a, input_sig)

    error = np.abs(output_lp3 - output_tf)
    max_error = np.max(error)

    plt.plot(input_sig, color="black", alpha=0.05, lw=1, label="input")
    plt.plot(output_lp3, color="blue", alpha=0.5, lw=1, label="lp3")
    plt.plot(output_chamberlin, color="orange", alpha=0.5, lw=1, label="chamberlin")
    # plt.plot(output_tf, color="red", alpha=0.5, lw=1, label="TF")
    # plt.plot(output_lp3 - output_tf)
    plt.grid()
    plt.legend()
    plt.title("Comparison of LP3 and Chambelin SVF")
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.show()


def plot_chamberlin_lp3_response(
    B_fixed: float, A_values: list, num_points: int = 1000
):
    w = np.geomspace(1 / 50000, np.pi, num_points)
    f_n = w / (2 * np.pi)
    z = np.exp(1j * w)

    plt.figure(figsize=(8, 4))

    colors = plt.cm.viridis(np.linspace(0, 0.85, len(A_values)))
    for idx, A in enumerate(A_values):
        # if np.isclose(B_fixed, 1.0):
        #     raise ValueError("B cannot be exactly 1.0.")

        numerator = (A / (1.0 - B_fixed)) * z * (z - B_fixed)
        denominator = z**2 + (A - B_fixed - 1.0) * z + B_fixed
        H = numerator / denominator

        magnitude_db = 20 * np.log10(np.abs(H) + np.finfo(np.float64).eps)
        plt.plot(f_n, magnitude_db, color=colors[idx], label=f"A = {A:.2e}")
        plt.axvline(A, color=colors[idx], alpha=0.66, lw=1, ls="--")

    plt.axvline(1 / 6, color="black", ls="--", lw=1, label=r"$f_n=1/6$")
    plt.axhline(20 * np.log10(1 / 2), color="black", ls="-.", lw=1, label=r"-6 dB")

    plt.title(f"ChamberlinLp3 Amplitude Response (B = {B_fixed:.3f})")
    plt.xlabel(r"Normalized Frequency [rad/2π]")
    plt.ylabel("Magnitude [dB]")
    plt.grid(which="both", color="#f8f8f8")
    plt.ylim([-20, 20])
    plt.xscale("log")
    plt.legend(loc="upper left")
    plt.tight_layout()
    # plt.savefig("chamberlin_lp3_amplitude_response.svg")
    plt.show()


if __name__ == "__main__":
    # plot_equivalence()
    plot_chamberlin_lp3_response(
        1 / np.sqrt(2),
        ema_alpha(np.geomspace(1e-5, 0.1, 9)),
        2048,
    )
