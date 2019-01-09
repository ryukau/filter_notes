import matplotlib.pyplot as pyplot
import numpy
import scipy.signal
import soundfile
from burgers import Burgers1D

samplerate = 44100
duration = 0.4
frequency = 60
amp = 0.8
dc = 0.6

phase = numpy.linspace(
    0,
    2 * numpy.pi * frequency * duration,
    int(samplerate * duration),
)
input_signal = numpy.sin(phase)

# 入力信号を整形。
amp = dc * amp if dc < 0.5 else (1 - dc) * amp
signal = (dc - amp) + amp * (input_signal + 1)

input_signal = numpy.copy(signal)

# シミュレーション。
burgers = Burgers1D(257)
read_index = burgers.wave.shape[1] - 2
for i, y in enumerate(signal):
    burgers.pick(0, y)
    burgers.step()
    signal[i] = burgers.state()[read_index]

raw_signal = numpy.copy(signal)

# ポップノイズを防ぐために信号の開始からしきい値を超えるまでの値を置き換え。
threshold = numpy.median(signal)
i = 0
while signal[i] < threshold:
    i += 1
signal[:i] = threshold
signal -= threshold

# 直流の除去。
hp_sos = scipy.signal.butter(4, 2 * 20 / samplerate, btype="highpass", output="sos")
signal = scipy.signal.sosfilt(hp_sos, signal)

time = [i / samplerate for i in range(len(signal))]

pyplot.title("Input Signal")
pyplot.xlabel("Time [s]")
pyplot.ylabel("Amplitude")
pyplot.ylim([0, 1])
pyplot.grid(color="#eeeeee")
pyplot.plot(time, input_signal, label="Signal")
pyplot.axhline(y=dc, lw=1, color="black", alpha=0.5, label="DC")
pyplot.axhline(y=dc - amp, lw=1, ls="--", color="black", alpha=0.5, label="Amp Low")
pyplot.axhline(y=dc + amp, lw=1, ls="--", color="black", alpha=0.5, label="Amp High")
pyplot.legend()
pyplot.show()

pyplot.title("Raw Output Signal")
pyplot.xlabel("Time [s]")
pyplot.ylabel("Amplitude")
pyplot.ylim([0, 1])
pyplot.grid(color="#eeeeee")
pyplot.plot(time, raw_signal, label="Signal")
pyplot.axhline(y=dc, lw=1, color="black", alpha=0.5, label="DC")
pyplot.axhline(y=dc - amp, lw=1, ls="--", color="black", alpha=0.5, label="Amp Low")
pyplot.axhline(y=dc + amp, lw=1, ls="--", color="black", alpha=0.5, label="Amp High")
pyplot.legend()
pyplot.show()

pyplot.title("Processed Output Signal")
pyplot.xlabel("Time [s]")
pyplot.ylabel("Amplitude")
pyplot.ylim([-1, 1])
pyplot.grid(color="#eeeeee")
pyplot.plot(time, signal, label="Signal")
pyplot.legend()
pyplot.show()
