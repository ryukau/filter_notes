import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import soundfile

data, fs = soundfile.read("snd/output.wav")

plt.plot(data)
plt.grid()
# plt.legend()
plt.show()
