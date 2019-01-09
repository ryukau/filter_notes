"""
- Range class
    https://stackoverflow.com/questions/12116685/how-can-i-require-my-python-scripts-argument-to-be-a-float-between-0-0-1-0-usin

- threshold_int function
    https://stackoverflow.com/questions/18700634/python-argparse-integer-condition-12
"""

import argparse
import functools
import multiprocessing
import numpy
import scipy.signal
import soundfile
import subprocess
from pathlib import Path

class Burgers1D:
    """
    dx = 1, dt = 1 で固定した inviscid Burgers' equation 。
    """

    def __init__(self, length):
        """length はシミュレーションする波の配列の長さ。"""
        self.wave = numpy.zeros((2, length))
        self.last = self.wave.shape[1] - 1

    def step(self, pick_y, read_index):
        self.wave = numpy.roll(self.wave, 1, axis=0)
        wave_m = self.wave[1]
        wave_l = numpy.roll(wave_m, 1)
        wave_r = numpy.roll(wave_m, -1)
        u_star_l = numpy.where(
            wave_m >= wave_r,
            numpy.where((wave_m + wave_r) * 0.5 > 0, wave_m, wave_r),
            numpy.where(
                wave_m > 0,
                wave_m,
                numpy.where(wave_r < 0, wave_r, 0),
            ),
        )
        u_star_r = numpy.where(
            wave_l >= wave_m,
            numpy.where((wave_l + wave_m) * 0.5 > 0, wave_l, wave_m),
            numpy.where(
                wave_l > 0,
                wave_l,
                numpy.where(wave_m < 0, wave_m, 0),
            ),
        )
        self.wave[0] = wave_m - (u_star_l * u_star_l - u_star_r * u_star_r) * 0.5
        self.wave[0][0] = pick_y
        self.wave[0][self.last] = 0
        return self.wave[0][read_index]

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def __eq__(self, other):
        return self.start <= other <= self.end

    def __repr__(self):
        return f"{self.start}-{self.end}"

def threshold_int(x, threshold):
    x = int(x)
    if x < threshold:
        raise argparse.ArgumentTypeError(
            f"Value must be greater or equal to {threshold}.")
    return x

def job_burgers(signal, length, read_index):
    output = numpy.empty_like(signal)
    burgers = Burgers1D(length)
    for i, pick_y in enumerate(signal):
        output[i] = burgers.step(pick_y, read_index)
    return output

def apply_burgers(src_sig, amp, dc, length, oversampling, raw=False):
    """
    src_sig.shape = (channel, frame)
    """
    amp = dc * amp if dc < 0.5 else (1 - dc) * amp
    src_sig = (dc - amp) + amp * (src_sig + 1)

    read_index = length - 2

    job_burgers_args = [(signal, length, read_index) for signal in src_sig]
    with multiprocessing.Pool() as pool:
        dest_sig = pool.starmap(job_burgers, job_burgers_args)
    dest_sig = numpy.array(dest_sig)

    if not raw:
        for sig in dest_sig:
            threshold = numpy.median(sig)
            i = 0
            while sig[i] < threshold:
                i += 1
            sig[:i] = threshold
            sig -= threshold

    return numpy.array(dest_sig)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Modify sound with inviscid Burgers' equation.")
    parser.add_argument("source", type=str, help="Input wav file path.")
    parser.add_argument(
        "-a",
        "--amp",
        type=float,
        default=1.0,
        choices=[Range(0.0, 1.0)],
        help="Amplitude of signal. Default to 1.0.")
    parser.add_argument(
        "-d",
        "--dc",
        type=float,
        default=0.5,
        choices=[Range(0.0, 1.0)],
        help="Direct current or offset from 0. Default to 0.5.")
    parser.add_argument(
        "-l",
        "--length",
        type=functools.partial(threshold_int, threshold=4),
        default=64,
        help="Length of wave or distance from audio source. Default to 64.")
    parser.add_argument(
        "-s",
        "--oversampling",
        type=functools.partial(threshold_int, threshold=1),
        default=4,
        help="Oversampling. Default to 4.")
    parser.add_argument(
        "--raw",
        action="store_true",
        help="If specified, output signal won't be normalized."
    )
    parser.add_argument(
        "-o",
        "--output",
        metavar="dest",
        type=str,
        help="Output file path. If not specified, add 'burgers_' as prefix to " \
             "input filename. For example, the output of 'path/to/some.wav' " \
             "will be 'path/to/burgers_some.wav'."
    )
    args = parser.parse_args()

    source = Path(str(args.source))
    dest = args.output
    if dest is None:
        dest = source.parent / Path("burgers_" + source.name)
    temp_path = source.parent / Path("🚧temp🚧.wav")
    amp = args.amp
    dc = args.dc
    length = args.length + 1  # +1 because boundary of wave is always 0.
    oversampling = args.oversampling
    raw = args.raw

    # sox was faster than python-samplerate. It might depends on environment.

    source_fs = subprocess.check_output(
        ['soxi', '-r', str(source)],
        encoding="utf-8",
    ).rstrip()
    source_fs = int(source_fs)
    temp_fs = source_fs * oversampling
    process = subprocess.run(
        ["sox", "--norm", source, temp_path, "rate", "-m",
         str(temp_fs)])
    if process.returncode != 0:
        exit(1)

    src_sig, fs = soundfile.read(str(temp_path), always_2d=True)
    dest_sig = apply_burgers(src_sig.T, amp, dc, length, oversampling, raw)
    soundfile.write(str(temp_path), dest_sig.T, temp_fs)

    command = ["sox", temp_path, dest, "rate", "-v", str(source_fs)]
    if not raw:
        command += ["highpass", "20", "gain", "-n"]
    process = subprocess.run(command)
    if process.returncode != 0:
        exit(1)

    temp_path.unlink()
