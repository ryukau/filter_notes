import numpy as np
import matplotlib.pyplot as pyplot

def envelope1(a, b, time):
    return (1 - np.exp(-a * time)) * np.exp(-b * time)

def envelope2(a, b, time):
    return (1 - a**time) * b**time

attack = 1.0
decay = 2.0
eps = 1e-5

a = -np.log(eps) / attack
b = -np.log(eps) / decay
peakTime = np.log(a / b + 1) / a
gain = 1 / envelope1(a, b, peakTime)

samplerate = 48000
duration = 1
time = np.linspace(0, duration, duration * samplerate)

atk = np.power(eps, 1 / attack)
dec = np.power(eps, 1 / decay)

out1 = envelope1(a, b, time)
# out2 = envelope2(atk, dec, time)

pyplot.figure(figsize=(6, 3))
pyplot.title(f"ExpAD: A={attack:3.1f}, B={decay:3.1f}")
pyplot.axvline(peakTime, color="black", alpha=0.5, lw=1, ls=":", label="t_p")
pyplot.plot(time, out1, label="Raw")
pyplot.plot(time, out1 * gain, label="Normalized")
# pyplot.plot(time, out2, label="out2")
pyplot.xlabel("Time [s]")
pyplot.ylabel("Amplitude")
pyplot.legend()
pyplot.grid()
pyplot.tight_layout()
pyplot.show()
