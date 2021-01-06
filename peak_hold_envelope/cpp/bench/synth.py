import numpy as np
import scipy.signal as signal
import soundfile

def nextTime(rng, rate):
    """
    Poisson process.
    https://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/

    If it happens once per N [someunit], then rate is set to 1/N.
    """
    return -np.log(1.0 - rng.uniform(0, 1)) / rate

def pulseNoise(rng, rate, length):
    out = np.zeros(length)
    time = 0
    while time < length:
        time += nextTime(rng, rate)

        t1 = int(time)
        t2 = t1 + 1
        amp = rng.uniform(0, 1)
        if t1 < length:
            out[t1] = amp * (t2 - time)
        if t2 < length:
            out[t2] = amp * (time - t1)
    return out

samplerate = 48000
rng = np.random.default_rng()
nFrames = 10 * samplerate
# nFrames = 1024

## Spiky noise.
sig = pulseNoise(rng, 1 / 8, nFrames)

## Worst case signal for peak detection.
# sig = rng.uniform(0, 1, nFrames)
# sig[::2] = 0

## Monotonically rising/falling
# sig = rng.uniform(0, 1, nFrames)
# sig = np.sort(sig)
# # sig = sig[::-1]  # Uncomment to render falling signal.

soundfile.write("snd/input.wav", sig, samplerate, subtype="FLOAT")
