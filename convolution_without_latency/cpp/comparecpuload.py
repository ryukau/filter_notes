import numpy as np
import matplotlib.pyplot as plt
import soundfile

length = 2**15
timeImmediate, fs = soundfile.read("snd/time_ImmediateConvolver.wav")
timeSplit, fs = soundfile.read("snd/time_SplitConvolver16_11.wav")

fig, ax = plt.subplots(2, 1)
ax[0].set_title("Concentrated Load")
ax[0].plot(timeImmediate[:length], color="black")
# ax[0].plot(timeImmediate[:length], color="blue", alpha=0.5, lw=1)
# ax[0].plot(timeSplit[:length], color="red", alpha=0.25, lw=1)
ax[0].set_ylabel("CPU Time Spent [ms]")

ax[1].set_title("Distributed Load")
ax[1].plot(timeSplit[:length], color="black")
ax[1].set_ylabel("CPU Time Spent [ms]")
ax[1].set_xlabel("Audio Time [sample]")
# ax[1].set_ylim((-0.002, 0.082))

for axis in ax:
    axis.set_ylim((-0.02, 1.02))
    axis.grid(color="#f8f8f8")

fig.set_size_inches((6, 6))
plt.tight_layout()
plt.show()
