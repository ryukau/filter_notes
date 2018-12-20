import functools
import multiprocessing

from pitch import *

def generate_sin(duration, samplerate, frequency):
    length = int(duration * samplerate)
    phase = numpy.linspace(0, 2 * numpy.pi * frequency * duration, length)
    return numpy.sin(phase)

def generate_sin_with_noise(duration, samplerate, frequency, noise_ratio):
    length = int(duration * samplerate)
    phase = numpy.linspace(0, 2 * numpy.pi * frequency * duration, length)
    signal = numpy.sin(phase) + noise_ratio * numpy.random.uniform(-1, 1)
    return signal / (noise_ratio + 1)

def generate_sin_am(duration, samplerate, car_freq, mod_freq):
    length = int(duration * samplerate)
    car_phase = numpy.linspace(0, 2 * numpy.pi * car_freq * duration, length)
    mod_phase = numpy.linspace(0, 2 * numpy.pi * mod_freq * duration, length)
    return numpy.sin(car_phase) * numpy.sin(mod_phase)

def generate_sin_fm(duration, samplerate, car_freq, mod_freq):
    length = int(duration * samplerate)
    mod_phase = numpy.linspace(0, 2 * numpy.pi * mod_freq * duration, length)
    car_phase = numpy.full(length, 2 * numpy.pi * car_freq * duration / length)
    phase = (car_phase + numpy.sin(mod_phase)).cumsum()
    return numpy.sin(phase)

class SignalData():
    def __init__(self, signal, samplerate, winlen, winstep, frequency, misc):
        self.signal = signal
        self.samplerate = samplerate
        self.winlen = winlen
        self.winstep = winstep
        self.frequency = frequency
        self.misc = misc

        self.pitch = {}
        for pitch_func in pitch_functions:
            self.pitch[pitch_func.__name__] = pitch_frame(
                signal, samplerate, winlen, winstep, pitch_func)

def test_sin_wave(samplerate, duration, winlen, winstep, freq_low, freq_high,
                  num):
    frequency = numpy.geomspace(freq_low, freq_high, num)
    signals = [generate_sin(duration, samplerate, freq) for freq in frequency]
    result = [
        SignalData(sig, samplerate, winlen, winstep, freq, None)
        for sig, freq in zip(signals, frequency)
    ]
    numpy.save("data/sin_wave.npy", result)

def job_sin_with_noise(samplerate, duration, winlen, winstep, frequency,
                       noise_ratio):
    return SignalData(
        generate_sin_with_noise(duration, samplerate, frequency, noise_ratio),
        samplerate,
        winlen,
        winstep,
        frequency,
        {"noise_ratio": noise_ratio},
    )

def test_sin_with_noise(samplerate, duration, winlen, winstep, freq_low,
                        freq_high, ratio_low, ratio_high, num):
    frequency = numpy.geomspace(freq_low, freq_high, num)
    noise_ratio = numpy.geomspace(ratio_low, ratio_high, num)
    result = []
    for freq in frequency:
        args = [(samplerate, duration, winlen, winstep, freq, ratio)
                for ratio in noise_ratio]
        with multiprocessing.Pool() as pool:
            result.append(pool.starmap(job_sin_with_noise, args))
    numpy.save(f"data/sin_with_noise.npy", result)

def job_sin_modulator(samplerate, duration, winlen, winstep, car_freq, mod_freq,
                      generation_func):
    return SignalData(
        generation_func(duration, samplerate, car_freq, mod_freq),
        samplerate,
        winlen,
        winstep,
        car_freq,
        {"mod_freq": mod_freq},
    )

def test_sin_am(samplerate, duration, winlen, winstep, car_freq_low,
                car_freq_high, mod_freq_low, mod_freq_high, num):
    car_freqs = numpy.geomspace(car_freq_low, car_freq_high, num)
    mod_freqs = numpy.geomspace(mod_freq_low, mod_freq_high, num)
    result = []
    for car_freq in car_freqs:
        args = [(
            samplerate,
            duration,
            winlen,
            winstep,
            car_freq,
            mod_freq,
            generate_sin_am,
        ) for mod_freq in mod_freqs]
        with multiprocessing.Pool() as pool:
            result.append(pool.starmap(job_sin_modulator, args))
    numpy.save(f"data/sin_am.npy", result)

def test_sin_fm(samplerate, duration, winlen, winstep, car_freq_low,
                car_freq_high, mod_freq_low, mod_freq_high, num):
    car_freqs = numpy.geomspace(car_freq_low, car_freq_high, num)
    mod_freqs = numpy.geomspace(mod_freq_low, mod_freq_high, num)
    result = []
    for car_freq in car_freqs:
        args = [(
            samplerate,
            duration,
            winlen,
            winstep,
            car_freq,
            mod_freq,
            generate_sin_fm,
        ) for mod_freq in mod_freqs]
        with multiprocessing.Pool() as pool:
            result.append(pool.starmap(job_sin_modulator, args))
    numpy.save(f"data/sin_fm.npy", result)

if __name__ == "__main__":
    samplerate = 16000
    duration = 0.2
    num = 32
    winlen = 0.1
    winstep = 0.01

    freq_low = 10
    freq_high = 8000

    # test_sin_wave(samplerate, duration, winlen, winstep, freq_low, freq_high,
    #               num,)

    # ノイズの平均をとるため duration を長くする。
    ratio_low = 0.01
    ratio_high = 100
    # test_sin_with_noise(samplerate, duration * 4, winlen, winstep, freq_low,
    #                     freq_high, ratio_low, ratio_high, num)

    mod_freq_low = 0.1
    mod_freq_high = 800
    # test_sin_am(samplerate, duration * 4, winlen, winstep, freq_low, freq_high,
    #             mod_freq_low, mod_freq_high, num)
    # test_sin_fm(samplerate, duration * 4, winlen, winstep, freq_low, freq_high,
    #             mod_freq_low, mod_freq_high, num)
