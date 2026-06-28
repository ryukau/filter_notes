import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class ComplexLowpass1:
    def __init__(self):
        self.x1 = 0
        self.y1 = 0

    def process(self, x0, cutoffNormalized, R):
        """
        cutoffNormalized in [0, 0.5).
        R in [0, 1].
        """
        cut = np.exp(2j * np.pi * cutoffNormalized)
        res = np.exp(cut.imag * np.log(R))

        a = cut * res
        b = (1 - a) / 2

        self.y1 = b * (x0 + self.x1) + a * self.y1
        self.x1 = x0
        return self.y1.real


def plot_response(R=0, cutoffs=np.geomspace(0.4 / 2**10, 0.4, 11)):
    N = 2**16
    eps = np.finfo(float).eps

    freqs = np.fft.rfftfreq(N)

    plt.figure(figsize=(8, 5))

    colors = cm.viridis(np.linspace(0, 0.85, len(cutoffs)))
    for i, fc in enumerate(cutoffs):
        flt = ComplexLowpass1()

        ir_input = np.zeros(N)
        ir_input[0] = 1.0
        ir_output = np.zeros(N)

        for n in range(N):
            ir_output[n] = flt.process(ir_input[n], fc, R)

        magnitude = np.abs(np.fft.rfft(ir_output))
        mag_db = 20 * np.log10(magnitude + eps)

        plt.plot(freqs, mag_db, label=f"Cutoff: {fc:.5f}", color=colors[i])
        plt.axvline(fc, color=colors[i], alpha=0.66, lw=1, ls="--")

    plt.title(f"1-pole Complex Lowpass Response (R={R:.4f})")
    plt.xlabel("Normalized Frequency [rad/2π]")
    plt.ylabel("Amplitude [dB]")

    plt.grid(True, color="#f0f0f0")
    plt.ylim(-20, 20)
    plt.xlim(2e-5, 0.5)
    plt.xscale("log")

    plt.legend(loc="upper left", frameon=True)
    plt.tight_layout()
    plt.savefig("img/amplitude_response.svg")
    plt.show()


def plot_ir():
    fc = 0.4
    R = 0.243
    N = 2**16

    svf = ComplexLowpass1()

    ir_input = np.zeros(N)
    ir_input[0] = 1.0
    ir_output = np.zeros(N)

    for n in range(N):
        ir_output[n] = svf.process(ir_input[n], fc, R)

    plt.figure(figsize=(10, 6))
    plt.plot(ir_output, color="black")
    plt.title(f"1-pole Complex Impulse Response (R={R})")
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.grid(True, color="#f0f0f0")
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    # plot_ir()
    plot_response(R=np.exp(-np.sqrt(2)))
    # plot_response(R=0.9)
