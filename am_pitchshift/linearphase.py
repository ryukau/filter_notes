"""
Proof of concept of low quality but linear phase frequency shifter.

Algorithm:

0. Prepare 2 buffers with same length.
1. Fill a buffer.
2. Apply frequency shift by using Hilbert transform (FFT is used here).
3. Crossfade into another buffer.
4. When the playback of a buffer is finished, go back to 1.

Combination of shifting amount and buffer length affects the quality of frequency shift. It seems like the quality becomes okay when following equation is satisfied.

```
bufferLength = sampleRateHz / shiftHz.
```
"""

import matplotlib.pyplot as plt
import numpy
import scipy.signal as signal
import soundfile
import sys

from pathlib import Path
from numpy.polynomial import polynomial


def add_delay(sos):
    return numpy.vstack((sos, [0, 1, 0, 1, 0, 0]))


def frequency_shift_mod(samplerate, analytic_signal, ratio):
    norm = numpy.abs(analytic_signal)
    theta = numpy.angle(analytic_signal)
    time = numpy.linspace(0, len(analytic_signal) / samplerate, len(analytic_signal))
    return norm * numpy.cos(ratio * theta + time)


def frequency_shift(samplerate, analytic_signal, shift_hz):
    norm = numpy.abs(analytic_signal)
    theta = numpy.angle(analytic_signal)
    time = numpy.linspace(0, len(analytic_signal) / samplerate, len(analytic_signal))
    return norm * numpy.cos(theta + 2 * numpy.pi * shift_hz * time)


def naive(samplerate, sig, shift_hz=1000):
    return frequency_shift(samplerate, signal.hilbert(sig), shift_hz)


def linear_phase_frequency_shift(samplerate, source, shift_hz, frameLength=512):
    frameLength += frameLength % 2  # Must be even.
    half = frameLength // 2
    length = frameLength * (len(source) // frameLength + 1) + half

    sig = numpy.zeros(length)
    sig[half : half + len(source)] = source

    window = numpy.hstack([numpy.linspace(0, 1, half), numpy.linspace(1, 0, half)])
    # window = signal.get_window("blackman", frameLength)
    dest = numpy.zeros(length)

    start = 0
    while start < len(sig) - half:
        end = start + frameLength

        frame = sig[start:end]
        dest[start:end] += window * naive(samplerate, frame, shift_hz)

        start += half
    return dest


def test_frequency_shift(samplerate, source, shift_hz=200):
    frameLength = 330

    out_lin = linear_phase_frequency_shift(samplerate, source, shift_hz, frameLength)
    soundfile.write("snd/linearphase.wav", out_lin, samplerate)

    out_ref = naive(samplerate, source, shift_hz)
    soundfile.write("snd/linearphase_reference.wav", out_ref, samplerate)


def test_mix(samplerate, source=None):
    """
    When `source is None`, the amplitude response should be all close to 0 dB.
    """
    frameLength = 512

    if source is None:
        source = numpy.zeros(2**16)
        source[0] = 1

    out_lin = linear_phase_frequency_shift(samplerate, source, 0, frameLength)

    half = frameLength // 2
    out_lin[half : half + len(source)] += source  # Change to `-=` to check difference.
    out_lin /= 2

    freq, response = signal.freqz(out_lin, worN=2048, fs=48000)

    gain = 20 * numpy.log10(numpy.maximum(numpy.abs(response), 1e-7))

    fig, ax = plt.subplots(2, 1)
    fig.set_size_inches((6, 8))

    ax[0].set_title("Amplitude Response")
    ax[0].plot(freq, gain)
    ax[0].set_xlabel("Frequency [Hz]")
    ax[0].set_ylabel("Amplitude [dB]")
    # ax[0].set_ylim([-6, 6])
    ax[0].set_xscale("log")

    ax[1].set_title("Signal")
    ax[1].plot(out_lin)
    ax[1].set_xlabel("Time [sample]")
    ax[1].set_ylabel("Amplitude")

    for axis in ax:
        axis.grid()
    plt.tight_layout()
    plt.show()

    soundfile.write("snd/linearphase_mix.wav", out_lin, samplerate)


data, samplerate = soundfile.read("snd/yey.wav", always_2d=True)
source = data.T[0]
print(samplerate)
test_frequency_shift(samplerate, source, 267)
test_mix(samplerate, None)
