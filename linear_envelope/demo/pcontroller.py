import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

def getKp(ω_c=np.pi / 16):
    kp = (np.sqrt(2) * np.exp(2 * 1j * ω_c) +
          (1 - np.sqrt(2)) * np.exp(1j * ω_c) - 1) / (2 * np.exp(2 * 1j * ω_c) - 1)

    # y = 1 - np.cos(ω_c)
    # kp = -y + np.sqrt((y + 2) * y)

    return kp

def getResponse(kp):
    b = np.array([kp, 0])
    a = np.array([1, kp - 1])
    return signal.freqz(b, a, 2**16)

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

def plotResponse():
    nKp = 16
    kpArray = np.linspace(1 / 100, 1, nKp)
    response = [getResponse(kp) for kp in kpArray]
    plotAmplitude(kpArray, response)
    plotPhase(kpArray, response)

def plotOmegaToKp():
    cmap = plt.get_cmap("plasma")
    plt.figure(figsize=(8, 4))
    omega_c = np.geomspace(20 / 48000, 20000 / 48000, 8)
    for idx, ω_c in enumerate(omega_c):
        kp = getKp(ω_c)
        response = getResponse(kp)
        color = cmap(idx / len(omega_c))
        plt.plot(
            response[0],
            20 * np.log10(abs(response[1])),
            color=color,
            lw=1,
            label=f"ω_c={ω_c:3f}")
        plt.axvline(ω_c, color=color, alpha=0.5, lw=1)
    plt.title("Frequency Responce and ω_c")
    plt.ylabel("Amplitude [dB]")
    plt.xlabel("Frequency [rad/sample]")
    plt.xscale("log")
    plt.ylim([-4, 1])
    plt.grid()
    plt.legend(loc=1)
    plt.tight_layout()
    plt.show()

plotResponse()
plotOmegaToKp()
