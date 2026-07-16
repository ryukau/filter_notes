import scipy.signal as signal
import numpy as np
import matplotlib.pyplot as plt
import math
import sys


class Hardclip:
    @staticmethod
    def f0(x):
        if x <= -1.0:
            return -1.0
        if x >= 1.0:
            return 1.0
        return x

    @staticmethod
    def f1(x):
        z = abs(x)
        if z < 1.0:
            return 0.5 * x * x
        return z - 0.5

    @staticmethod
    def f2(x):
        if abs(x) < 1.0:
            return (x * x * x) / 6.0
        s = 1.0 if x >= 0 else -1.0
        return (s * x - 1) * x / 2 + s / 6

    @staticmethod
    def f3(x):
        x2 = x * x
        if x2 < 1.0:
            return (x2 * x2) / 24.0
        s = 1.0 if x >= 0 else -1.0
        return (((s * 4 * x - 6) * x + s * 4) * x - 1) / 24

    @staticmethod
    def f4(x):
        z = abs(x)
        if z < 1.0:
            return (x**5) / 120.0
        s = 1.0 if x >= 0 else -1.0
        return s * (((((z - 2.0) * z + 2.0) * z - 1.0) * z / 24.0) + 1.0 / 120.0)

    @staticmethod
    def g1(x):
        if abs(x) >= 1.0:
            return 0.0
        return 1.0

    @staticmethod
    def g2(x):
        return 0.0


class NaiveSaturator:
    """
    A helper class for comparison.
    """

    def __init__(self, nonlinear_func):
        self.fn = nonlinear_func

    def reset(self):
        pass

    def process_array(self, input_sig):
        output_sig = np.zeros_like(input_sig)
        for i, val in enumerate(input_sig):
            output_sig[i] = self.fn.f0(val)
        return output_sig


class Adaa_DividedDifference:
    def __init__(self, order, funcs):
        """funcs.fN are N-th antiderivatives. funcs.gN are N-th derivatives."""
        self.order = order
        self.funcs = funcs
        self.x = [0.0] * order
        self.f = [0.0] * order

        e = math.sqrt(sys.float_info.epsilon)
        self.threshold = (1920 * e) ** (1 / 5)

        self.reset()

    def _call_f(self, n, val):
        return getattr(self.funcs, f"f{n}")(val)

    def _call_correction(self, n, val):
        return self._call_f(n, val) if n >= 0 else getattr(self.funcs, f"g{-n}")(val)

    def reset(self, input_val=0.0):
        self.x = [input_val] * self.order
        self.f = [self._call_f(self.order - i, input_val) for i in range(self.order)]

    def process(self, input_val):
        curr = self._call_f(self.order, input_val)
        next_f = [0.0] * self.order

        for step in range(self.order):
            target = self.order - 1 - step
            next_f[step] = curr
            d = input_val - self.x[step]

            if abs(d) < self.threshold:
                mid = (input_val + self.x[step]) * 0.5
                curr = self._call_f(target, mid) + (
                    (d * d) / 24.0
                ) * self._call_correction(target - 2, mid)
            else:
                curr = (step + 1) * (curr - self.f[step]) / d

        self.f = next_f
        self.x = [input_val] + self.x[:-1]
        return curr

    def process_array(self, input_sig):
        output_sig = np.zeros_like(input_sig)
        for i, val in enumerate(input_sig):
            output_sig[i] = self.process(val)
        return output_sig


def run_comparison_test(use_window=False):
    sr = 48000
    duration = 0.5
    # gain = 8.0
    gain = 1e3

    freq_spec = 1234.56
    freq_wave = 2.0

    cmap = plt.get_cmap("hot")
    cases = [
        ("Naive", NaiveSaturator(Hardclip), 0.2),
        ("Adaa DD 1", Adaa_DividedDifference(1, Hardclip), 0.3),
        ("Adaa DD 2", Adaa_DividedDifference(2, Hardclip), 0.3),
        ("Adaa DD 3", Adaa_DividedDifference(3, Hardclip), 0.3),
        ("Adaa DD 4", Adaa_DividedDifference(4, Hardclip), 0.3),
    ]

    # --- Processing ---
    t = np.linspace(0, duration, int(sr * duration), endpoint=False)
    sig_spec_in = gain * np.sin(2 * np.pi * freq_spec * t)
    sig_wave_in = gain * np.sin(2 * np.pi * freq_wave * t)

    results_spec = []
    results_wave = []
    for name, prc, _ in cases:
        out_s = prc.process_array(sig_spec_in)
        results_spec.append(out_s)

        prc.reset()
        out_w = prc.process_array(sig_wave_in)
        results_wave.append(out_w)

    # --- Plotting ---
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))

    ax1.set_title(f"Frequency Domain (Freq {freq_spec}Hz, Gain {gain})")
    ax1.set_xlabel("Frequency (Hz)")
    ax1.set_ylabel("Magnitude (dBFS)")
    ax1.set_ylim(-120, 10)
    ax1.grid(True, which="both", alpha=0.3)
    window = signal.get_window("hann", len(sig_spec_in), fftbins=True)
    for i, (name, _, alpha) in enumerate(cases):
        sig = results_spec[i]
        if use_window:
            sig *= window
        fft = np.fft.fft(sig, norm="forward")
        mag = 20 * np.log10(np.abs(fft[: len(sig) // 2]) + 1e-16)
        freqs = np.linspace(0, sr / 2, len(mag))

        color = cmap(i / len(cases))
        ax1.plot(freqs, mag, label=name, color=color, alpha=alpha, linewidth=1.2)
    ax1.legend(loc="upper right")

    ax2.set_title(f"Time Domain (Freq {freq_wave}Hz, Gain {gain})")
    ax2.set_xlabel("Time (samples)")
    ax2.set_ylabel("Amplitude")
    ax2.set_ylim(-2, 2)
    ax2.grid(True, which="both", alpha=0.3)
    for i, (name, _, alpha) in enumerate(cases):
        color = cmap(i / len(cases))
        ax2.plot(results_wave[i], label=name, color=color, alpha=alpha, linewidth=1.2)
    ax2.legend(loc="upper right")

    plt.tight_layout()
    plt.savefig("img/adaa.svg")
    plt.show()


if __name__ == "__main__":
    run_comparison_test()
