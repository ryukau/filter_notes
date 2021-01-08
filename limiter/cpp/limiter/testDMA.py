"""
Test DoubleAverageFilter.
"""

import matplotlib.pyplot as plt
import numpy as np
import soundfile

sig, _ = soundfile.read("snd/input.wav")
reference, _ = soundfile.read("snd/smoothed.wav")
cpp, _ = soundfile.read("output/input_out.wav")

np.testing.assert_allclose(cpp, reference, rtol=1e-5)

delayed = np.hstack((np.zeros(32), sig[:-32]))
plt.plot(delayed, label="Input")
plt.plot(reference, label="Ref.")
plt.plot(cpp, label="C++")
plt.grid()
plt.legend()
plt.savefig("img/DoubleAverageFilterTest.svg")
