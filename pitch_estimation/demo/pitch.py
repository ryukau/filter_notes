"""
Based on the code written by sevagh.

https://github.com/sevagh/pitch-detection/blob/master/src/yin.cpp
"""

import numpy
import python_speech_features
from pyfftw.interfaces.numpy_fft import fft, ifft

YIN_THRESHOLD = 0.3

MPM_K = 0.5  # Type I NSD では後半に大きなピークができるので K を小さめに設定する。

### Common functions ###
def autocorrelation_type1(sig):
    len_sig = len(sig)
    spec = fft(sig)
    return ifft(spec * spec.conj()).real

def autocorrelation_type2(sig):
    len_sig = len(sig)
    sig = numpy.pad(sig, (0, len_sig), "constant")
    spec = fft(sig)
    return ifft(spec * spec.conj()).real[:len_sig]

def parabolic_interpolation(array, x):
    x_result = None
    if x < 1:
        x_result = x if array[x] <= array[x + 1] else x + 1
    elif x >= len(array) - 1:
        x_result = x if array[x] <= array[x - 1] else x - 1
    else:
        denom = array[x + 1] + array[x - 1] - 2 * array[x]
        delta = array[x - 1] - array[x + 1]
        if denom == 0:
            return x
        return x + delta / (2 * denom)
    return x_result

### YIN ###
def difference_type1(sig):
    autocorr = autocorrelation_type1(sig)
    return autocorr[0] - autocorr

def difference_type2(sig):
    autocorr = autocorrelation_type2(sig)
    energy = (sig * sig)[::-1].cumsum()[::-1]
    return energy[0] + energy - 2 * autocorr

def cumulative_mean_normalized_difference(diff):
    diff[0] = 1
    sum_value = 0
    for tau in range(1, len(diff)):
        sum_value += diff[tau]
        diff[tau] /= sum_value / tau
    return diff

def absolute_threshold(diff, threshold=YIN_THRESHOLD):
    tau = 2
    while tau < len(diff):
        if diff[tau] < threshold:
            while tau + 1 < len(diff) and diff[tau + 1] < diff[tau]:
                tau += 1
            break
        tau += 1
    return None if tau == len(diff) or diff[tau] >= threshold else tau

def invert_nsd(nsd):
    tau = 0
    while nsd[tau] > 0:
        nsd[tau] = 0
        tau += 1
        if tau >= len(nsd):
            return None
    return -nsd

def yin_cmnd_type1(sig, samplerate):
    diff = difference_type1(sig)
    cmnd = cumulative_mean_normalized_difference(diff)
    tau = absolute_threshold(cmnd)
    if tau is None:
        return numpy.nan
    return samplerate / parabolic_interpolation(cmnd, tau)

def yin_cmnd_type2(sig, samplerate):
    diff = difference_type2(sig)
    cmnd = cumulative_mean_normalized_difference(diff)
    tau = absolute_threshold(cmnd)
    if tau is None:
        return numpy.nan
    return samplerate / parabolic_interpolation(cmnd, tau)

def yin_nsd_type1(sig, samplerate):
    nsd = normalized_square_difference_type1(sig)
    nsd = invert_nsd(nsd)
    if nsd is None:
        return numpy.nan
    tau = absolute_threshold(nsd, 0)
    if tau is None:
        return numpy.nan
    return samplerate / tau

def yin_nsd_type2(sig, samplerate):
    nsd = normalized_square_difference_type2(sig)
    nsd = invert_nsd(nsd)
    if nsd is None:
        return numpy.nan
    tau = absolute_threshold(nsd, 0)
    if tau is None:
        return numpy.nan
    return samplerate / tau

### MPM ###
def normalized_square_difference_type1(sig):
    corr = autocorrelation_type1(sig)
    return corr / corr[0] if corr[0] != 0 else None

def normalized_square_difference_type2(sig):
    corr = autocorrelation_type2(sig)
    cumsum = (sig * sig)[::-1].cumsum()[::-1]
    cumsum[cumsum < 1] = 1  # 発散を防ぐ。
    return corr / (corr[0] + cumsum)

def estimate_period(diff):
    start = 0
    while diff[start] > 0:
        start += 1
        if start >= len(diff):
            return None

    threshold = MPM_K * numpy.max(diff[start:])
    isNegative = True
    max_index = None
    for i in range(start, len(diff)):
        if isNegative:
            if diff[i] < 0:
                continue
            max_index = i
            isNegative = False
        if diff[i] < 0:
            isNegative = True
            if diff[max_index] >= threshold:
                return max_index
        if diff[i] > diff[max_index]:
            max_index = i
    return None

def mpm_cmnd_type1(sig, samplerate):
    diff = difference_type1(sig)
    cmnd = cumulative_mean_normalized_difference(diff)
    cmnd = numpy.max(cmnd) / 2 - cmnd
    index = estimate_period(cmnd)
    if index is None:
        return numpy.nan
    return samplerate / parabolic_interpolation(cmnd, index)

def mpm_cmnd_type2(sig, samplerate):
    diff = difference_type2(sig)
    cmnd = cumulative_mean_normalized_difference(diff)
    cmnd = numpy.max(cmnd) / 2 - cmnd
    index = estimate_period(cmnd)
    if index is None:
        return numpy.nan
    return samplerate / parabolic_interpolation(cmnd, index)

def mpm_nsd_type1(sig, samplerate):
    nsd = normalized_square_difference_type1(sig)
    index = estimate_period(nsd)
    if index is None:
        return numpy.nan
    return samplerate / parabolic_interpolation(nsd, index)

def mpm_nsd_type2(sig, samplerate):
    nsd = normalized_square_difference_type2(sig)
    index = estimate_period(nsd)
    if index is None:
        return numpy.nan
    return samplerate / parabolic_interpolation(nsd, index)

def framesig(data, frame_len, frame_step):
    """
    検証用。
    """
    n_frame = int(len(data) / frame_step)
    frame = numpy.zeros((n_frame, frame_len))
    for i in range(n_frame):
        start = i * frame_step
        end = start + frame_len
        part = data[start:end]
        frame[i][:len(part)] = part[:]
    return frame

def pitch_frame(data, samplerate, winlen, winstep, pitch_func=yin_cmnd_type2):
    frame = python_speech_features.sigproc.framesig(
        data,
        frame_len=int(samplerate * winlen),
        frame_step=int(samplerate * winstep),
    )
    return numpy.array([pitch_func(sig, samplerate) for sig in frame])

pitch_functions = [
    yin_cmnd_type1,
    yin_cmnd_type2,
    yin_nsd_type1,
    yin_nsd_type2,
    mpm_cmnd_type1,
    mpm_cmnd_type2,
    mpm_nsd_type1,
    mpm_nsd_type2,
]
