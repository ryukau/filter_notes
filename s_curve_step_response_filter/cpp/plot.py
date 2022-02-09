import soundfile
import numpy as np
import matplotlib.pyplot as plt

data, sampleRate = soundfile.read("snd/bessel.wav")

plt.plot(data, label="Output")
plt.axvline(len(data) // 2, color="black", alpha=0.5, lw=1, ls="--", label="Target Time")
plt.grid()
plt.legend()
plt.show()
