import argparse
import functools
import matplotlib.pyplot as pyplot
import numpy
import python_speech_features
import soundfile
import scipy.signal
from multiprocessing import Pool
from pathlib import Path
from pyfftw.interfaces.numpy_fft import fft, ifft

def trim_feature(feature, n_frame):
    feature_step = feature.shape[0]
    if feature_step < n_frame:
        zeros = numpy.zeros((n_frame - feature_step, feature.shape[1]))
        feature = numpy.concatenate((feature, zeros))
    elif feature_step > n_frame:
        feature = feature[0:n_frame]
    return feature

def get_envelope(data, samplerate, winstep, n_frame):
    data_abs = numpy.abs(data)
    envelope = numpy.zeros(n_frame)
    index = None
    step = int(winstep * samplerate)
    start = 0
    end = step
    for frame in range(n_frame):
        if end >= len(data_abs):
            index = frame
            break
        envelope[frame] = numpy.max(data_abs[start:end])
        start = end
        end += step
    if index is not None:
        envelope[index] = numpy.max(data_abs[start:])
        index += 1
    return envelope

def trim_data(data, n_frame):
    if len(data) < n_frame:
        return numpy.pad(data, (0, n_frame - len(data)), "constant")
    elif len(data) > n_frame:
        return data[:n_frame]
    return data

def autocorrelation_type2(sig):
    len_sig = len(sig)
    sig = numpy.pad(sig, (0, len_sig), "constant")
    spec = fft(sig)
    return ifft(spec * spec.conj()).real[:len_sig]

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

def normalized_square_difference_type2(sig):
    corr = autocorrelation_type2(sig)
    cumsum = (sig * sig)[::-1].cumsum()[::-1]
    cumsum[cumsum < 1] = 1  # 発散を防ぐ。
    return corr / (corr[0] + cumsum)

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

MPM_K = 0.5

def get_key_maxima(diff):
    start_index = 0
    while start_index < len(diff) and diff[start_index] > 0:
        start_index += 1
    isNegative = True
    max_index = 0
    key_maxima = []
    for i in range(start_index, len(diff)):
        if isNegative:
            if diff[i] < 0:
                continue
            max_index = i
            isNegative = False
        if diff[i] < 0:
            isNegative = True
            key_maxima.append((max_index, diff[max_index]))
        if diff[i] > diff[max_index]:
            max_index = i
    return key_maxima  # [(index, value), ...]

def hz_to_cent(frequency):
    return 1200 * numpy.log2(frequency / 440)

def mpm_nsd_type2(sig, samplerate, n_pitch):
    nsd = normalized_square_difference_type2(sig)
    key_maxima = get_key_maxima(nsd)

    key_maxima = [
        hz_to_cent(samplerate / parabolic_interpolation(nsd, index))
        for index, _ in sorted(key_maxima, key=lambda x: x[1], reverse=True)
    ]
    if len(key_maxima) < n_pitch:
        key_maxima += [-1000] * (n_pitch - len(key_maxima))
    else:
        key_maxima = key_maxima[0:n_pitch]
    return sorted(key_maxima)

def pitch_frame(data, samplerate, winlen, winstep, n_pitch, pitch_func):
    frame = python_speech_features.sigproc.framesig(
        data,
        frame_len=int(samplerate * winlen),
        frame_step=int(samplerate * winstep),
    )
    return numpy.array([pitch_func(sig, samplerate, n_pitch) for sig in frame])

def extract_feature(data, samplerate, winstep, nfft, n_frame):
    winlen = nfft / samplerate
    numcep = 26
    n_pitch = 18

    mfcc = trim_feature(
        python_speech_features.mfcc(
            data,
            samplerate,
            winlen=winlen,
            winstep=winstep,
            numcep=numcep,
            nfilt=numcep * 2,
            nfft=nfft,
            lowfreq=0,
            highfreq=20000,
            preemph=0.0,
            ceplifter=0,
        ),
        n_frame,
    )
    mfcc_max = numpy.max(mfcc)
    if mfcc_max > 0:
        mfcc /= mfcc_max

    delta = python_speech_features.base.delta(mfcc, 1)
    envelope = get_envelope(data, samplerate, winstep, n_frame)

    pitch = pitch_frame(data, samplerate, winlen, winstep, n_pitch,
                        mpm_nsd_type2)
    pitch = numpy.ravel(pitch)
    pitch[numpy.isnan(pitch)] = 0
    pitch = numpy.concatenate((
        pitch,
        numpy.zeros(n_frame * n_pitch - len(pitch)),
    ))

    return (
        envelope,
        numpy.concatenate((numpy.ravel(mfcc), numpy.ravel(delta))),
        pitch,
    )

class Sound:
    def __init__(self, filepath, envelope, mfccdelta, pitch):
        self.filepath = filepath
        self.envelope = envelope
        self.mfccdelta = mfccdelta
        self.pitch = pitch

def get_feature(directory_path, nfft=1024, n_frame=100):
    winstep = 0.01  # Frame length in second.

    filepath = [path for path in directory_path.glob("*.wav")]
    sound = [soundfile.read(str(path)) for path in filepath]

    extract = functools.partial(
        extract_feature, winstep=winstep, nfft=nfft, n_frame=n_frame)

    with Pool() as pool:
        features = list(pool.starmap(extract, sound))

    return [Sound(path, *data) for path, data in zip(filepath, features)]

def plot_envelope(filepath, envelope):
    for env, path in zip(envelope, filepath):
        pyplot.plot(env, label=path.stem)
    pyplot.legend()
    pyplot.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract features from wav file.")
    parser.add_argument(
        "src_name",
        metavar="src_name",
        type=str,
        help="source directory name. snd/src_name")
    args = parser.parse_args()

    source_name = args.src_name
    source_directory = Path("snd") / Path(source_name)
    if not source_directory.exists():
        print(f"{str(source_directory)} does not exist")
        exit()

    Path("data").mkdir(exist_ok=True)

    sounds = get_feature(source_directory)

    # for snd in sounds:
    #     print(len(snd.pitch))

    numpy.save(f"data/{source_name}.npy", sounds)
