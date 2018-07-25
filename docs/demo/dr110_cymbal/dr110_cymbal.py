import argparse
import datetime
import multiprocessing
import numpy as np
import pathlib
import soundfile
import scipy.signal as signal


def normalize(data, peak=1.0):
    return peak * data / np.abs(np.max(data))


def pulse(time, frequency):
    noise = np.power(10, -(np.random.uniform(0, 5, len(time)) + 5))
    return signal.square(2 * np.pi * frequency * (time + noise))


def noise(size):
    return np.random.uniform(-1, 1, size)


def randomAmp(low, high):
    lowLog = np.log(low)
    diffLog = np.log(high) - lowLog
    return np.exp(lowLog + diffLog * np.random.random())


def cymbalSource(time, pulseFreq=[317, 465, 820, 1150]):
    sig = 0
    for freq in pulseFreq:
        sig += randomAmp(0.999, 1) * pulse(time, freq)
    return sig + noise(len(time)) / 3.3


def applyContinuousFilter(samplerate, source, system):
    num, den, dt = signal.cont2discrete(system, 1.0 / samplerate)
    return signal.lfilter(num[0], den, source)


def highMetalFilter(samplerate, source):
    return applyContinuousFilter(samplerate, source, (
        [9.471e+10, 2.5625e+16, 0],
        [3437973.0, 1.816815e+11, 8.03e+15, 3.125e+20],
    ))


def lowMetalFilter(samplerate, source):
    return applyContinuousFilter(samplerate, source, (
        [2.68345e+11, 3.5234375e+16, 0],
        [30108309.0, 5.227075e+11, 1.56671875e+16, 1.953125e+20],
    ))


def rcHighpass(samplerate, source, r, c):
    rc = r * c
    num, den, dt = signal.cont2discrete(([rc, 0], [rc, 1]), 1.0 / samplerate)
    return signal.lfilter(num[0], den, source)


def rcLowpass(samplerate, source, r, c):
    num, den, dt = signal.cont2discrete(([1], [r * c, 1]), 1.0 / samplerate)
    return signal.lfilter(num[0], den, source)


def tHighpass(samplerate, source, r, c1, c2):
    num, den, dt = signal.cont2discrete(
        ([c1 * c2 / (c1 + c2), 0, 0], [1, 1, 1 / r / (c1 + c2)]),
        1.0 / samplerate,
    )
    return signal.lfilter(num[0], den, source)


def envelope(time, T, E, threshold=0.4):
    a = np.log(threshold / E) / T * 1000
    return np.exp(a * time)


highMetalEnvTime = {
    "OH": 700,
    "CH": 80,
    "CY": 60,
}


def metalMix(samplerate, duration, metalType="CY", detune=0,
             envThreshold=0.4):
    time = np.linspace(0, duration, int(samplerate * duration), endpoint=False)

    source = cymbalSource(time, [
        317 * (1 + np.random.uniform(-detune, detune)),
        465 * (1 + np.random.uniform(-detune, detune)),
        820 * (1 + np.random.uniform(-detune, detune)),
        1150 * (1 + np.random.uniform(-detune, detune)),
    ])

    high_metal = highMetalFilter(samplerate, source)
    high_metal = rcHighpass(samplerate, high_metal, 1e6, 47e-9)
    high_metal = rcHighpass(samplerate, high_metal, 5600, 470e-12)
    high_metal = rcHighpass(samplerate, high_metal, 100e3, 10e-10)
    high_metal = rcHighpass(samplerate, high_metal, 5600, 22e-10)
    env567 = envelope(time, highMetalEnvTime.get(metalType, "CY"), 6,
                      envThreshold)
    high_metal_mix = high_metal * env567
    if (metalType == "CY"):
        env8 = envelope(time, 900, 6, envThreshold)
        high_metal_mix += high_metal * env8 / 10
    high_metal = normalize(high_metal_mix)

    low_metal = 0
    if (metalType == "CY"):
        low_metal = lowMetalFilter(samplerate, source)
        low_metal = rcHighpass(samplerate, low_metal, 1e6, 47e-9)
        low_metal = rcHighpass(samplerate, low_metal, 68000, 470e-12)
        low_metal = rcHighpass(samplerate, low_metal, 10000, 470e-12)
        env9 = envelope(time, 1400, 2.7, envThreshold)
        low_metal = low_metal * env9 / 10
        low_metal = normalize(low_metal) * randomAmp(0.4, 0.5)

    return normalize(low_metal + high_metal)


def smoothOut(sig):
    """
    Smoothing out the artifacts at the end of sample that is caused by
    np.signal.resample.
    """
    size = int(len(sig) * 0.005)
    offset = len(sig) - size
    factor = np.pi / (size - 1) / 2
    for i in range(0, size):
        index = i + offset
        sig[index] *= np.cos(i * factor)
    return sig


def renderMetal(filepath, samplerate, oversampling, duration, metalType,
                detune, envThreshold):
    sig = metalMix(samplerate, duration, metalType, detune, envThreshold)
    sig = signal.resample(sig, int(len(sig) / oversampling))
    if oversampling != 1:
        sig = smoothOut(sig)
    soundfile.write(filepath, sig, int(samplerate / oversampling))


if __name__ == "__main__":
    output_dir = "snd"
    num_render_per_type = 8
    oversampling = 16
    detune = 0.01
    samplerate = 44100 * oversampling

    time_stamp = datetime.datetime.now().strftime(
        "DR-110_Cymbals_%Y-%m-%d-%H-%M-%S-%f")
    output_dir = pathlib.Path(output_dir) / pathlib.Path(time_stamp)
    output_dir.mkdir(parents=True, exist_ok=True)

    args = []
    durationDict = {"OH": 1, "CH": 0.2, "CY": 4}
    for metalType in highMetalEnvTime:
        for i in range(num_render_per_type):
            filename = "{}_{:04d}_{}.wav".format(metalType, i, time_stamp)
            args.append((
                str(pathlib.Path(output_dir) / pathlib.Path(filename)),
                samplerate,
                oversampling,
                durationDict[metalType],
                metalType,
                detune,
                0.4 * np.exp(-i / 2),
            ))

    pool = multiprocessing.Pool()
    pool.starmap_async(renderMetal, args)
    pool.close()
    pool.join()
