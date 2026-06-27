import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class ChamberlinSVF:
    def __init__(self):
        self.lp = 0
        self.bp = 0

    def process_fast(self, x0, cutoffNormalized, Q):
        cut = np.clip(cutoffNormalized, 0, 0.1731886233119285)
        f = 2 * np.sin(np.pi * cut)

        q = np.clip(1 / Q, 0, 2 - f)

        hp = x0 - self.lp - self.bp * q
        self.bp += f * hp
        self.lp += f * self.bp
        return self.lp

    def process_full_range(self, x0, cutoffNormalized, Q):
        cut = np.clip(cutoffNormalized, 0, 0.4997)
        f = 2 * np.sin(np.pi * cut)

        q = np.clip(1 / Q, 0, ((2 - f) * (2 + f)) / (2 * f))

        hp = x0 - self.lp - self.bp * q
        self.bp += f * hp
        self.lp += f * self.bp
        return self.lp


def plot_min_Q():
    cut = np.geomspace(1e-5, 0.5 - 1e-5, 2**16)
    f = 2 * np.sin(np.pi * cut)

    min_Q_exact = 2 * f / ((2 - f) * (2 + f))

    x = 2 - f
    min_Q_fast_1 = 1 / (2 - f)
    min_Q_fast_2 = 1 / (x * (1 + x * 1 / 4))
    min_Q_fast_3 = 1 / (x * (1 + x * (1 / 4 + x * 1 / 8)))

    colors = plt.cm.plasma(np.linspace(0, 1, 3, endpoint=False))

    plt.figure(figsize=(8, 5))
    plt.plot(
        cut,
        min_Q_exact,
        label=r"Exact  $\frac{2f}{4 - f^2}$",
        color="black",
        alpha=0.25,
        lw=1,
    )
    plt.plot(
        cut,
        min_Q_fast_1,
        label=r"$1\ /\ \left[ (2 - f) \right]$",
        color=colors[0],
        alpha=0.75,
        lw=1,
    )
    plt.plot(
        cut,
        min_Q_fast_2,
        label=r"$1\ /\ \left[ (2 - f) + \frac{(2 - f)^2}{2^2} \right]$",
        color=colors[1],
        alpha=0.75,
        lw=1,
    )
    plt.plot(
        cut,
        min_Q_fast_3,
        label=r"$1\ /\ \left[ (2 - f) + \frac{(2 - f)^2}{2^2} + \frac{(2 - f)^3}{2^3} \right]$",
        color=colors[2],
        alpha=0.75,
        lw=1,
    )
    plt.axhline(
        1 / np.sqrt(2),
        color="red",
        ls="--",
        lw=1,
        alpha=0.5,
        label=r"Max. flat Q  $\frac{1}{\sqrt{2}}$",
    )
    plt.axvline(
        np.asin((np.sqrt(6) - np.sqrt(2)) / 2) / np.pi,
        color="blue",
        ls="--",
        lw=1,
        alpha=0.5,
        label=r"Max. flat bound ~ 0.173",
    )
    plt.title(f"Chamberlin SVF, Cutoff vs Minimum Q")
    plt.xlabel("Normalized Cutoff [rad/2π]")
    plt.ylabel("Min. Q")
    # plt.xscale("log")
    plt.yscale("log")
    plt.ylim([1e-2, 1e4])
    plt.grid(which="both", color="#f0f0f0")
    plt.legend()
    plt.tight_layout()
    # plt.savefig("img/cutoff_vs_min_Q.svg")
    plt.show()


def plot_chamberlin_response(Q=1.0, cutoffs=np.geomspace(0.4 / 2**10, 0.4, 11)):
    N = 2**16
    eps = np.finfo(float).eps

    freqs = np.fft.rfftfreq(N)

    plt.figure(figsize=(8, 5))

    colors = cm.viridis(np.linspace(0, 0.85, len(cutoffs)))
    for i, fc in enumerate(cutoffs):
        svf = ChamberlinSVF()

        ir_input = np.zeros(N)
        ir_input[0] = 1.0
        ir_output = np.zeros(N)

        for n in range(N):
            ir_output[n] = svf.process_fast(ir_input[n], fc, Q)
            # ir_output[n] = svf.process_full_range(ir_input[n], fc, Q)

        magnitude = np.abs(np.fft.rfft(ir_output))
        mag_db = 20 * np.log10(magnitude + eps)

        plt.plot(freqs, mag_db, label=f"Cutoff: {fc:.5f}", color=colors[i])
        plt.axvline(fc, color=colors[i], alpha=0.66, lw=1, ls="--")

    monotonic_ub = np.asin(np.sqrt(np.sqrt(6) - np.sqrt(2)) / 2) / np.pi
    plt.axvline(monotonic_ub, color="black", ls="--", lw=1, label="Max. flat bound")
    plt.axhline(20 * np.log10(1 / 4), color="black", ls="-.", lw=1, label=r"-12 dB")

    plt.title(f"Chamberlin SVF Lowpass Response (Q={Q:.4f})")
    plt.xlabel("Normalized Frequency [rad/2π]")
    plt.ylabel("Amplitude [dB]")

    plt.grid(True, color="#f0f0f0")
    plt.ylim(-20, 20)
    plt.xlim(2e-5, 0.5)
    plt.xscale("log")

    plt.legend(loc="upper left", frameon=True)
    plt.tight_layout()
    # plt.savefig("img/amplitude_response.svg")
    plt.show()


def plot_chamberlin_ir():
    fc = 0.4
    Q = 16
    N = 2**16

    svf = ChamberlinSVF()

    ir_input = np.zeros(N)
    ir_input[0] = 1.0
    ir_output = np.zeros(N)

    for n in range(N):
        ir_output[n] = svf.process_full_range(ir_input[n], fc, Q)

    plt.figure(figsize=(10, 6))
    plt.plot(ir_output, color="black")
    plt.title(f"Chamberlin SVF Impulse Response (Q={Q})")
    plt.xlabel("Time [sample]")
    plt.ylabel("Amplitude")
    plt.grid(True, color="#f0f0f0")
    plt.legend(frameon=True)
    plt.tight_layout()
    plt.show()


def save_parabola_svg(filename):
    fig, ax = plt.subplots()

    offset = 0.3
    x0 = np.linspace(-1, 1, 2048) + offset
    x1 = np.linspace(-1, 1, 2048)
    x2 = np.linspace(-1, 1, 2048) - offset
    ax.plot(x0, -((x0 - offset) ** 2), color="#2288ff", ls=":")
    ax.plot(x1, -(x1**2), color="black")
    ax.plot(x2, -((x2 + offset) ** 2), color="#ff6688", ls=":")

    ax.set_axis_off()
    plt.savefig(
        filename, format="svg", transparent=True, bbox_inches="tight", pad_inches=0
    )
    # plt.close(fig)
    plt.show()


def plot_P():
    def P(ω, f, Q):
        a1 = f * f + f / Q - 2
        a2 = 1 - f / Q
        return (
            4 * a2 * np.cos(ω) ** 2
            + (2 * a1 * a2 + 2 * a1) * np.cos(ω)
            + a2**2
            - 2 * a2
            + a1**2
            + 1
        )

    f_array = np.linspace(0, 1, 11)
    colors = plt.cm.plasma(np.linspace(0, 1, len(f_array), endpoint=False))
    for i, f in enumerate(f_array):
        x = np.linspace(0, np.pi, 1024)
        y = P(x, f, 1 / np.sqrt(2))
        plt.plot(x, y, color=colors[i], label=f"{f:.1f}")
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    # plot_min_Q()
    # plot_chamberlin_ir()
    plot_chamberlin_response(Q=1 / np.sqrt(2))
    # save_parabola_svg("parabola.svg")
    # plot_P()
