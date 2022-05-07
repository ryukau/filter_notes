import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import soundfile

out, fs = soundfile.read("output/oracleengine_out.wav")
peak, fs = soundfile.read("output/oracleengine_peak.wav")
smoothed, fs = soundfile.read("output/oracleengine_smoothed.wav")
delayed, fs = soundfile.read("output/oracleengine_delayed.wav")

absed = np.abs(out)
maxIndex = np.argmax(absed)
print(maxIndex, out[maxIndex])

absed = absed[maxIndex - 200:maxIndex + 200]
peak = peak[maxIndex - 200:maxIndex + 200]

# fir = signal.get_window("bartlett", int(0.002 * fs))
# fir /= np.sum(fir)
# fir /= np.sum(fir)
# inv = np.where(peak > 1, 1 / peak, 1)
# smoothed = np.convolve(inv, fir)[:-len(fir)]
smoothed = smoothed[maxIndex - 200:maxIndex + 200]
delayed = delayed[maxIndex - 200:maxIndex + 200]

# absed = absed[:1024]
# peak = peak[:1024]
# smoothed = smoothed[:1024]
# delayed = delayed[:1024]

# plt.plot(absed, lw=1, alpha=0.5, label="absed")
# plt.plot(peak, lw=1, alpha=0.5, label="peak")
plt.plot(1 / smoothed, lw=1, alpha=0.5, label="smoothed")
plt.plot(np.abs(delayed), lw=1, alpha=0.5, label="delayed")
plt.grid(color="#f0f0f0")
plt.legend()
plt.show()
