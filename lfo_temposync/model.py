import numpy as np
import matplotlib.pyplot as plt

fs = 48000
tempo = 120
sync0 = 0.25
sync1 = 1.0

duration = 1
nFrame = int(fs * duration)

p0 = 0.5
p1 = 0.0

nLfo = 0.1 * fs
mid = nLfo / 2

v0 = tempo / (fs * sync0 * 60)
v1 = tempo / (fs * sync1 * 60)

distance = p1 + v1 * nLfo - p0
k = np.ceil((v0 + v1) / 2 * mid - distance)
height = (distance + k) / mid - (v0 + v1) / 2

velocity = np.interp(np.arange(nFrame), [0, mid, nLfo, nFrame], [v0, height, v1, v1])

phase1 = np.linspace(p1, p1 + v1 * nFrame, nFrame) % 1
phase2 = (p0 + velocity.cumsum()) % 1

time = np.linspace(0, duration, nFrame)
plt.figure(figsize=(6, 3), tight_layout=True)
plt.plot(time, phase1, label="Target (p1)")
plt.plot(time, phase2, label="Smoothed (p2)")
plt.xlabel("Time [s]")
plt.ylabel("Phase")
plt.title("Comparison of Target Phase and Smoothed Phase")
plt.grid()
plt.legend(loc=4)
plt.show()
# plt.savefig("img/PythonImplementationFixed.svg")
