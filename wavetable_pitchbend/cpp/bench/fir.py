import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
import json

def splitPolyPhase(fir, nPhase):
    padding = nPhase - len(fir) % nPhase
    fir = np.hstack((fir, np.zeros(padding)))

    polyphase = []
    for phase in range(nPhase):
        part = fir[phase::nPhase]
        polyphase.append(part.tolist())

    return polyphase

def getDiff(polyphase):
    poly = np.array(polyphase)
    rolled = np.roll(poly, -1, 0)
    rolled[-1] = np.roll(rolled[-1], 1)
    return rolled - poly

def formatFir(polyphase):
    text = ""
    for phase in polyphase:
        text += "{"
        for value in phase:
            text += f"Sample({value}),"
        text = text[:-1]
        text += "},\n"
    return text

def decimationFir(length, oversample, cutoff, stop=0.5, fs=48000, plot=True):
    bands = np.hstack((
        [0, cutoff],
        [stop, oversample / 2],
    ))
    desired = (1, 0)
    weight = (1, 100)
    fir = signal.remez(length, bands, desired, weight, fs=oversample, maxiter=1024)

    if plot:
        freq, h0 = signal.freqz(fir, worN=2048, fs=oversample * fs)
        plt.plot(freq, ampToDecibel(h0), label=f"{oversample}x")
        # plt.xscale("log")
        plt.grid(which="both")
        plt.legend()
        plt.show()

    return fir

def toCpp(fir, N, M):
    poly = splitPolyPhase(M * fir, M)[::-1]
    polyText = formatFir(poly)

    diff = getDiff(poly)
    diffText = formatFir(diff)

    text = """\
#pragma once

#include <array>

template<typename Sample> struct TableInterpFir {
"""
    text += f"constexpr static size_t bufferSize = {N};"
    text += f"constexpr static size_t nPhase = {M};"
    text += """
  constexpr static std::array<std::array<Sample, bufferSize>, nPhase> coefficient{{
"""
    text += polyText
    text += """\
  }};

  constexpr static std::array<std::array<Sample, bufferSize>, nPhase> diff{{
"""
    text += diffText
    text += """\
  }};
};
"""
    with open("fir.hpp", "w", encoding="utf-8") as fi:
        fi.write(text)

f_cp = 0.4  # cutoff
f_cs = 0.5
N = 32  # number of taps per phase.
M = 64  # number of phase.
fir = decimationFir(N * M - 1, M, f_cp, f_cs, fs=1, plot=False)

toCpp(fir, N, M)
