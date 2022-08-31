import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt
import soundfile
import json

def firlp2hp(lp):
    """
    `index` search is heuristic. Check if something is wrong.
    """
    hp = -lp
    index = np.argmax(np.abs(hp))
    hp[index] += 1
    return hp

def generateSin(frequency, duration, sampleRate=48000):
    length = int(sampleRate * duration)
    phase = np.linspace(0, 2 * np.pi * frequency * duration, length)
    return (np.sin(phase), sampleRate)

nTap = 2**15 - 1

# cpp, fs = soundfile.read("snd/ImmediateConvolver.wav")
# time, fs = soundfile.read("snd/time_ImmediateConvolver.wav")

cpp, fs = soundfile.read("snd/SplitConvolver16_11.wav")
time, fs = soundfile.read("snd/time_SplitConvolver16_11.wav")

add = soundfile.read("snd/FFTConvolutionFilter.wav")[0][nTap:]

coefficient = signal.firwin(nTap, 1000, window="nuttall", fs=48000)
coefficient = np.hstack([firlp2hp(coefficient), [0]])
# coefficient[:len(coefficient) // 2] = 0
# coefficient[len(coefficient) // 4 * 3:] = 0

# sigA, fs = generateSin(1000, 4)
# sigB, fs = generateSin(100, 4)
# sig = sigA + sigB
# sig = sigB

sig = np.hstack([np.ones(int(4 * fs // 2)), np.zeros(int(4 * fs // 2))])

scipy_convolve = signal.convolve(sig, coefficient)

plt.figure()
plt.title("Output")
# plt.plot(sig, alpha=0.5, label="input")
plt.plot(cpp, alpha=0.5, label="cpp")
# plt.plot(add, alpha=0.5, label="add")
plt.plot(scipy_convolve, alpha=0.5, label="SciPy")
plt.grid()
plt.legend()

plt.figure()
plt.title("Diff")
plt.plot(cpp - scipy_convolve[:len(cpp)], alpha=0.5, label="diff_cpp")
# plt.plot(add - scipy_convolve[:len(add)], alpha=0.5, label="diff_add")
plt.grid()
plt.legend()

plt.figure()
plt.title("Time")
plt.plot(time, alpha=0.5, label="time")
plt.xlabel("Time [sample]")
plt.ylabel("Time Spent [ms]")
plt.grid()
plt.legend()

plt.show()
