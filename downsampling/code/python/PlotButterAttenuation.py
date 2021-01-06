import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from itertools import cycle

def ampToDecibel(z):
    return 20 * np.log10(np.abs(z))

orderList = np.arange(4, 17)
oversampleList = [2, 4, 8, 16]
attenuationFreq = [0, 0.75, 0.875]

cutoff = 0.4

markerIter = cycle(["1", "2", "3", "4"])
fig, ax = plt.subplots(1, 2)
fig.set_size_inches(8, 4)
fig.set_tight_layout(True)
fig.suptitle(f"Butterworth Filter Aliasing Attenuation (f_c=0.4)")

cmap = plt.get_cmap("viridis")
for plotIdx, oversample in enumerate(oversampleList):
    attenuation = [[] for _ in range(len(attenuationFreq) - 1)]
    for order in orderList:
        sos = signal.butter(
            order,
            cutoff,
            "low",
            output="sos",
            fs=oversample,
        )
        freq, h0 = signal.sosfreqz(sos, worN=attenuationFreq, fs=oversample)
        attn = [int(x) for x in ampToDecibel(h0)][1:]
        for idx, value in enumerate(attn):
            attenuation[idx].append(value)

    color = cmap(plotIdx / len(oversampleList))
    marker = next(markerIter)
    ax[0].scatter(orderList,
                  attenuation[0],
                  color=color,
                  marker=marker,
                  label=f"{oversample}x")
    ax[1].scatter(orderList,
                  attenuation[1],
                  color=color,
                  marker=marker,
                  label=f"{oversample}x")

ax[0].set_title("Aliasing at f_s/4")
ax[1].set_title("Aliasing at f_s/8")
for axis in ax:
    axis.grid()
    axis.legend()
    axis.set_ylabel("Amplitude [dB]")
    axis.set_xlabel("Order")
    axis.set_ylim((-120, 5))

plt.show()
