import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sici

def blit(time):
    return np.sinc(time)

def blep(time):
    si, _ = sici(time * np.pi)
    return si / np.pi + 1 / 2

def blepResidual(time):
    step = np.where(time < 0, 0, np.where(time > 0, 1, 0.5))
    return blep(time) - step

def polyBlep4(t):
    if t < -2:
        return 0
    if t < -1:
        t += 2
        return t**4 / 24
    if t < 0:
        t += 1
        return -t**4 / 8 + t**3 / 6 + t**2 / 4 + t / 6 + 1 / 24
    if t < 1:
        return t**4 / 8 - t**3 / 3 + (2 * t) / 3 - 1 / 2
    if t < 2:
        t -= 1
        return -t**4 / 24 + t**3 / 6 - t**2 / 4 + t / 6 - 1 / 24
    return 0

def polyBlep6(t):
    if t < -3:
        return 0
    if t < -2:
        t += 3
        return t**6 / 720
    if t < -1:
        t += 2
        return -t**6 / 144 + t**5 / 120 + t**4 / 48 + t**3 / 36 + t**2 / 48 + t / 120 + 1 / 720
    if t < 0:
        t += 1
        return t**6 / 72 - t**5 / 30 - t**4 / 24 + t**3 / 18 + (5 * t**2) / 24 + (
            13 * t) / 60 + 29 / 360
    if t < 1:
        return -t**6 / 72 + t**5 / 20 - t**3 / 6 + (11 * t) / 20 - 1 / 2
    if t < 2:
        t -= 1
        return t**6 / 144 - t**5 / 30 + t**4 / 24 + t**3 / 18 - (5 * t**2) / 24 + (
            13 * t) / 60 - 29 / 360
    if t < 3:
        t -= 2
        return -t**6 / 720 + t**5 / 120 - t**4 / 48 + t**3 / 36 - t**2 / 48 + t / 120 - 1 / 720
    return 0

def polyBlep8(t):
    if t < -4:
        return 0
    if t < -3:
        t += 4
        return t**8 / 40320
    if t < -2:
        t += 3
        return -t**8 / 5760 + t**7 / 5040 + t**6 / 1440 + t**5 / 720 + t**4 / 576 + t**3 / 720 + t**2 / 1440 + t / 5040 + 1 / 40320
    if t < -1:
        t += 2
        return t**8 / 1920 - t**7 / 840 - t**6 / 360 + t**4 / 72 + t**3 / 30 + (
            7 * t**2) / 180 + t / 42 + 31 / 5040
    if t < 0:
        t += 1
        return -t**8 / 1152 + t**7 / 336 + t**6 / 288 - t**5 / 80 - (
            19 * t**4) / 576 + t**3 / 48 + (49 * t**2) / 288 + (397 *
                                                                t) / 1680 + 4541 / 40320
    if t < 1:
        return t**8 / 1152 - t**7 / 252 + t**5 / 45 - t**3 / 9 + (151 * t) / 315 - 1 / 2
    if t < 2:
        t -= 1
        return -t**8 / 1920 + t**7 / 336 - t**6 / 288 - t**5 / 80 + (
            19 * t**4) / 576 + t**3 / 48 - (49 * t**2) / 288 + (397 *
                                                                t) / 1680 - 4541 / 40320
    if t < 3:
        t -= 2
        return t**8 / 5760 - t**7 / 840 + t**6 / 360 - t**4 / 72 + t**3 / 30 - (
            7 * t**2) / 180 + t / 42 - 31 / 5040
    if t < 4:
        t -= 3
        return -t**8 / 40320 + t**7 / 5040 - t**6 / 1440 + t**5 / 720 - t**4 / 576 + t**3 / 720 - t**2 / 1440 + t / 5040 - 1 / 40320
    return 0

def polyBlepResidual(time, basisFunc):
    return np.array([basisFunc(t) for t in time])

nSample = 5

time = np.linspace(-nSample, nSample, 1024)

residual4 = polyBlepResidual(time, polyBlep4)
residual6 = polyBlepResidual(time, polyBlep6)
residual8 = polyBlepResidual(time, polyBlep8)

# plt.plot(time, blit(time))
# plt.plot(time, blep(time))
# plt.plot(time, blepResidual(time))

plt.title("PolyBLEP Residual")
plt.plot(time, residual4, lw=1, alpha=0.75, color="red", label="4 pt.")
plt.plot(time, residual6, lw=1, alpha=0.75, color="orange", label="6 pt.")
plt.plot(time, residual8, lw=1, alpha=0.75, color="blue", label="8 pt.")
plt.xlabel("Time [samples]")
plt.ylabel("Amplitude")
plt.legend()
plt.grid()
plt.show()
