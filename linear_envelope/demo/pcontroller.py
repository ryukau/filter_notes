import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

def getKp(omega_c=np.pi / 16):
    """正しい値にならない。"""
    # kp = -(np.sqrt(8 * np.exp(2j * omega_c) + (2**(5 / 2) - 8) * np.exp(1j * omega_c) - 2**
    #                (3 / 2) + 3) - np.sqrt(2) - 1) / (4 * np.exp(1j * omega_c) + 2**(3 / 2))
    kp = (np.sqrt(8 * np.exp(2j * omega_c) +
                  (2**(5 / 2) - 8) * np.exp(1j * omega_c) - 2**(3 / 2) + 3) + np.sqrt(2) +
          1) / (4 * np.exp(1j * omega_c) + 2**(3 / 2))
    # kp = (np.sqrt(8 * np.exp(2j * omega_c) +
    #               (-2**(5 / 2) - 8) * np.exp(1j * omega_c) + 2**(3 / 2) + 3) - np.sqrt(2) +
    #       1) / (4 * np.exp(1j * omega_c) - 2**(3 / 2))
    # kp = -(np.sqrt(8 * np.exp(2j * omega_c) + (-2**(5 / 2) - 8) * np.exp(1j * omega_c) + 2**
    #                (3 / 2) + 3) + np.sqrt(2) - 1) / (4 * np.exp(1j * omega_c) - 2**(3 / 2))
    return kp

def getResponse(kp):
    b = np.array([kp, 0])
    a = np.array([1, kp - 1])
    return signal.freqz(b, a, 8192)

def plotAmplitude(kpArray, response):
    cmap = plt.get_cmap("plasma")
    for idx, (w, h) in enumerate(response):
        plt.plot(
            w,
            20 * np.log10(abs(h)),
            color=cmap(idx / len(response)),
            label=f"k_p={kpArray[idx]:.3f}")
    plt.ylabel("Amplitude [dB]")
    plt.xlabel("Frequency [rad/sample]")
    plt.xscale("log")
    plt.grid()
    plt.legend()
    plt.show()

def plotPhase(kpArray, response):
    cmap = plt.get_cmap("plasma")
    for idx, (w, h) in enumerate(response):
        angles = np.unwrap(np.angle(h))
        plt.plot(
            w, angles, color=cmap(idx / len(response)), label=f"k_p={kpArray[idx]:.3f}")
    plt.ylabel("Angle (radians)")
    plt.xlabel("Frequency [rad/sample]")
    plt.xscale("log")
    plt.grid()
    plt.legend()
    plt.show()

nKp = 16
kpArray = np.linspace(1 / 100, 1, nKp)
response = [getResponse(kp) for kp in kpArray]
plotAmplitude(kpArray, response)
plotPhase(kpArray, response)
