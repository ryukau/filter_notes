import matplotlib.pyplot as pyplot
import numpy
import scipy.signal as signal
import soundfile
import sys

from pathlib import Path
from numpy.polynomial import polynomial

def plot_response(name, sos_real, sos_imag):
    _, h_real = signal.sosfreqz(sos_real, whole=True)
    w, h_imag = signal.sosfreqz(sos_imag, whole=True)

    h = h_real / h_imag
    # h = h_imag

    fig, ax1 = pyplot.subplots()
    ax1.set_title(f"Frequency response ({name})")
    ax1.plot(w, abs(h), "b", alpha=0.2)
    ax1.set_ylabel("Amplitude", color="b")
    ax1.set_xlabel("Frequency [rad/sample]")
    ax2 = ax1.twinx()
    angles = numpy.unwrap(numpy.angle(h))
    ax2.plot(w, angles, "g")
    ax2.set_ylabel("Angle (radians)", color="g")
    ax2.grid()
    ax2.axis("tight")

def add_delay(sos):
    return numpy.vstack((sos, [0, 1, 0, 1, 0, 0]))

def pitch_shift_mod(samplerate, analytic_signal, ratio):
    norm = numpy.abs(analytic_signal)
    theta = numpy.angle(analytic_signal)
    time = numpy.linspace(0, len(wav) / samplerate, len(wav))
    return norm * numpy.cos(ratio * theta + time)

def pitch_shift(samplerate, analytic_signal, shift_hz):
    norm = numpy.abs(analytic_signal)
    theta = numpy.angle(analytic_signal)
    time = numpy.linspace(0, len(wav) / samplerate, len(wav))
    return norm * numpy.cos(theta + shift_hz * time)

def naive(samplerate, sig, shift_hz=1000):
    return pitch_shift(samplerate, signal.hilbert(wav), shift_hz)

def niemitalo(samplerate, sig, shift_hz=1000):
    def section(a):
        a2 = a * a
        return [a2, 0, -1, 1, 0, -a2]

    sos_real = numpy.array([
        section(a) for a in [
            0.4021921162426,
            0.8561710882420,
            0.9722909545651,
            0.9952884791278,
        ]
    ])
    sos_imag = numpy.array([
        section(a) for a in [
            0.6923878000000,
            0.9360654322959,
            0.9882295226860,
            0.9987488452737,
        ]
    ])
    sos_imag = add_delay(sos_imag)

    # plot_response(sys._getframe().f_code.co_name, sos_real, sos_imag)

    real = signal.sosfilt(sos_real, sig)
    imag = signal.sosfilt(sos_imag, sig)
    analytic = 0.5 * (real + 1j * imag)
    return pitch_shift(samplerate, analytic, shift_hz)

def wasabi(samplerate, sig, shift_hz=1000):
    sos_real = numpy.array([
        [0.190696, 0, -1, 1, 0, -0.190696],
        [0.860735, 0, -1, 1, 0, -0.860735],
    ])
    sos_imag = numpy.array([[0.553100, 0, -1, 1, 0, -0.553100]])
    sos_imag = add_delay(sos_imag)

    # plot_response(sys._getframe().f_code.co_name, sos_real, sos_imag)

    real = signal.sosfilt(sos_real, sig)
    imag = signal.sosfilt(sos_imag, sig)
    analytic = 0.5 * (real + 1j * imag)
    return pitch_shift(samplerate, analytic, shift_hz)

def favreau(samplerate, sig, shift_hz=1000):
    def biquad(a1, a2):
        return [a2, a1, 1, 1, a1, a2]

    sos_real = numpy.array([
        biquad(0.02569, -0.260502),
        biquad(-1.8685, 0.870686),
    ])
    sos_imag = numpy.array([
        biquad(-1.94632, 0.94657),
        biquad(-0.83774, 0.06338),
    ])

    # plot_response(sys._getframe().f_code.co_name, sos_real, sos_imag)

    real = signal.sosfilt(sos_real, sig)
    imag = signal.sosfilt(sos_imag, sig)
    analytic = 0.5 * (real + 1j * imag)
    return pitch_shift(samplerate, analytic, shift_hz)

def wilkinson(samplerate, sig, shift_hz=1000):
    k = 9
    n = numpy.arange((k + 1) / 2)

    zero_real = numpy.exp(numpy.pi / 2**(2 * n))
    zero_real = numpy.append(zero_real, -zero_real)

    zero_imag = numpy.exp(numpy.pi / 2**(2 * n + 1))
    zero_imag = numpy.append(zero_imag, -zero_imag)

    sos_real = signal.zpk2sos(zero_real, 1 / zero_real, 1e-5)
    sos_imag = signal.zpk2sos(zero_imag, 1 / zero_imag, 1e-5)
    sos_imag = add_delay(sos_imag)

    # plot_response(sys._getframe().f_code.co_name, sos_real, sos_imag)

    real = signal.sosfilt(sos_real, sig)
    imag = signal.sosfilt(sos_imag, sig)
    analytic = real - 1j * imag

    sig = pitch_shift(samplerate, analytic, shift_hz)
    peak = numpy.max(numpy.abs(sig))
    return sig / peak if peak != 0 else sig

def allpass(samplerate, rc):
    rc *= 2 * numpy.pi
    num, den, dt = signal.cont2discrete(
        ([rc, -1], [rc, 1]),
        1 / samplerate,
        "gbt",
        0.5,
    )
    return num[0], den

def mcnulty(samplerate, sig, shift_hz=1000):
    rc_real = [
        allpass(samplerate, rc) for rc in [
            9.31e-06,
            4.2723e-05,
            0.0001836,
            0.00078146,
            0.003333,
            0.026055,
        ]
    ]
    rc_imag = [
        allpass(samplerate, rc) for rc in [
            2.6676e-06,
            2.08e-05,
            8.87e-05,
            0.00038064,
            0.0016049999999999999,
            0.007412,
        ]
    ]

    sos_real = [
        signal.tf2sos(
            polynomial.polymul(rc1[0], rc2[0]),
            polynomial.polymul(rc1[1], rc2[1]),
        )[0] for rc1, rc2 in zip(rc_real[::2], rc_real[1::2])
    ]
    sos_imag = [
        signal.tf2sos(
            polynomial.polymul(rc1[0], rc2[0]),
            polynomial.polymul(rc1[1], rc2[1]),
        )[0] for rc1, rc2 in zip(rc_imag[::2], rc_imag[1::2])
    ]

    # plot_response(sys._getframe().f_code.co_name, sos_real, sos_imag)

    real = signal.sosfilt(sos_real, sig)
    imag = signal.sosfilt(sos_imag, sig)

    analytic = 0.5 * (real - 1j * imag)
    return pitch_shift(samplerate, analytic, shift_hz)

def chuck(samplerate, sig, shift_hz=1000):
    rc_imag = [
        allpass(samplerate, rc)
        for rc in [5.49e-06, 4.75e-05, 2.37e-04, 1.27e-03]
    ]
    rc_real = [
        allpass(samplerate, rc)
        for rc in [2.00e-05, 1.07e-04, 5.36e-04, 4.64e-03]
    ]

    sos_imag = [
        signal.tf2sos(
            polynomial.polymul(rc1[0], rc2[0]),
            polynomial.polymul(rc1[1], rc2[1]),
        )[0] for rc1, rc2 in zip(rc_imag[::2], rc_imag[1::2])
    ]

    sos_real = [
        signal.tf2sos(
            polynomial.polymul(rc1[0], rc2[0]),
            polynomial.polymul(rc1[1], rc2[1]),
        )[0] for rc1, rc2 in zip(rc_real[::2], rc_real[1::2])
    ]

    # plot_response(sys._getframe().f_code.co_name, sos_real, sos_imag)

    real = signal.sosfilt(sos_real, sig)
    imag = signal.sosfilt(sos_imag, sig)

    analytic = 0.5 * (real - 1j * imag)
    return pitch_shift(samplerate, analytic, shift_hz)

# for wav in Path("snd").glob("*.wav"):
#     print(str(wav))

data, samplerate = soundfile.read("snd/yey.wav", always_2d=True)
wav = data.T[0]
shift_hz = 1000  # Hz

out = naive(samplerate, wav, shift_hz)
soundfile.write("snd/naive.wav", out, samplerate)

out = niemitalo(samplerate, wav, shift_hz)
soundfile.write("snd/niemitalo.wav", out, samplerate)

out = wasabi(samplerate, wav, shift_hz)
soundfile.write("snd/wasabi.wav", out, samplerate)

out = favreau(samplerate, wav, shift_hz)
soundfile.write("snd/favreau.wav", out, samplerate)

out = wilkinson(samplerate, wav, shift_hz)
soundfile.write("snd/wilkinson.wav", out, samplerate)

out = mcnulty(samplerate, wav, shift_hz)
soundfile.write("snd/mcnulty.wav", out, samplerate)

out = chuck(samplerate, wav, shift_hz)
soundfile.write("snd/chuck.wav", out, samplerate)

# pyplot.show()
