import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def blackmanHarris(length: int):
    isEven = 1 - length % 2

    ω = 2 * np.pi / float(length - isEven)
    φ = np.pi / 2
    k = 2 * np.cos(ω)
    u1 = np.sin(φ - ω)
    u2 = np.sin(φ - 2 * ω)

    window = np.zeros(length)
    for i in range(length):
        u0 = k * u1 - u2
        u2 = u1
        u1 = u0
        window[i] = 0.21747 + u0 * (-0.45325 + u0 * (0.28256 + u0 * -0.04672))
    return window


def testBlackmanHarris():
    length = 8

    ref = signal.get_window("blackmanharris", length, fftbins=False)
    fast = blackmanHarris(length)

    plt.plot(ref, label="ref.")
    plt.plot(fast, label="fast")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


def triangle(length: int):
    tri_N: int = length + 2
    inv_N: float = 1 / tri_N
    window = np.zeros(length)
    for i in range(length):
        window[i] = 1 - abs((2 * (i + 1) - tri_N) * inv_N)
    return window


def testTriangle():
    length = 32

    ref = signal.get_window("triang", length)
    fast = triangle(length)

    plt.plot(ref, label="ref.")
    plt.plot(fast, label="fast")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


def bartlett(length: int):
    tri_N: int = length
    inv_N: float = 1 / tri_N
    window = np.zeros(length)
    for i in range(length):
        window[i] = 1 - abs((2 * i - tri_N) * inv_N)
    return window


def triangle1(length: int):
    tri_N: int = length + 1
    inv_N: float = 1 / tri_N
    window = np.zeros(length)
    for i in range(length):
        window[i] = 1 - abs((2 * (i + 1) - tri_N) * inv_N)
    return window


def testBartlett():
    length = 31

    ref = signal.get_window("bartlett", length)
    fast = bartlett(length)

    plt.plot(ref, label="ref.")
    plt.plot(fast, label="fast")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    testBlackmanHarris()
    # testTriangle()
    # testBartlett()
