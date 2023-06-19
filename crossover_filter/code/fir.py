import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt


def kaiserExp(length, alpha):
    K = length // 2 + length % 2
    p = (np.arange(1, length + 1) - K) / (K)
    return np.exp(alpha * np.sqrt(1 - p * p)) / np.exp(alpha)


def kaiserCosh(length, alpha):
    K = length // 2 + length % 2
    p = (np.arange(1, length + 1) - K) / (K)
    return np.cosh(alpha * np.sqrt(1 - p * p)) / np.cosh(alpha)


def testKaiser():
    length = 15
    alpha = 10
    win_scipy = signal.get_window(("kaiser", np.pi * alpha), length + length % 2)

    # alpha = 1.3
    win_exp = kaiserExp(length, np.pi * alpha)
    win_cosh = kaiserCosh(length, np.pi * alpha)

    if length % 2 == 1:
        win_scipy = win_scipy[1:]
    else:
        win_scipy = win_scipy[1:]
        win_exp = win_exp[:-1]
        win_cosh = win_cosh[:-1]

    print(len(win_scipy), len(win_cosh))
    print(win_scipy[0], win_scipy[-1])
    print(win_cosh[0], win_cosh[-1])

    plt.plot(win_scipy, label="scipy")
    plt.plot(win_exp, label="exp")
    plt.plot(win_cosh, label="cosh")
    plt.grid()
    plt.legend()
    plt.show()


def testHann():
    length = 16
    win_scipy = signal.get_window("hann", length + length % 2)[::-1]

    K = length + length % 2
    n = np.arange(1, length + 1)
    sn = np.sin(np.pi * n / K)
    win_hann = sn * sn

    print(win_scipy[0], win_scipy[-1])

    plt.plot(win_scipy, label="scipy")
    plt.plot(win_hann, label="naive")
    plt.grid()
    plt.legend()
    plt.show()


def testLowpass():
    length = 15
    cutoff = 0.1
    lp = np.zeros(length)
    K = length // 2 + length % 2
    for idx in range(length):
        x = np.pi * (idx + 1 - K)
        if x != 0:
            lp[idx] = np.sin(cutoff * x * 2) / x
        else:
            lp[idx] = cutoff
    plt.plot(lp)
    plt.show()
