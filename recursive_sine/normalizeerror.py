import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special
import scipy.interpolate as interpolate
from numpy.polynomial import Polynomial


def approximatePolynomial(px: np.array, py: np.array):
    exponent = np.arange(len(px))
    A = np.array([x**exponent for x in px])
    return np.linalg.solve(A, py)


def derivative(degree, n):
    with np.errstate(divide="ignore"):
        return np.where(
            degree - n < 0, 0, special.factorial(degree) / special.factorial(degree - n)
        )


def approximatePolynomial2(px: np.array, py: np.array):
    degree = np.arange(len(px) + 1)
    A = np.array([x**degree for x in px])

    # first derivative
    s0 = np.arange(len(degree))
    s0[1] = 0
    A = np.vstack([A, s0])
    py = np.hstack([py, [0]])

    return np.linalg.solve(A, py)


def softclip(y):
    absed = np.abs(y)
    return np.where(absed <= 1, y, np.sign(y) * (np.tanh(absed - 1) + 1))


def generateRandomPolynomial(rng, nPoint, resolution=65536):
    px = rng.integers(1, resolution - 1, nPoint).astype(np.float64)
    px.sort()
    for i in range(len(px) - 1):
        if px[i] == px[i + 1]:
            offset = 1
            while i + offset < len(px) and px[i] == px[i + offset]:
                px[i + offset] += offset
                offset += 1
    px /= resolution

    py = rng.uniform(-1, 1, nPoint)

    px = np.hstack([[0], px, [1]])
    py = np.hstack([[0], py, [0]])
    poly = Polynomial(approximatePolynomial(px, py))
    return (poly, px, py)


def testRandomPolynomial():
    rng = np.random.default_rng(456)
    poly, px, py = generateRandomPolynomial(rng, 7)

    x = np.linspace(0, 1, 2**16)
    plt.scatter(px, py, zorder=2)

    y = poly(x)
    # y = softclip(y)
    plt.plot(x, y)
    plt.grid()
    plt.show()


rng = np.random.default_rng(156116)
nSample = 2 ** np.arange(3, 16)
nData = 1024
result = []
for seed in rng.integers(0, np.iinfo(np.uint64).max, nData, dtype=np.uint64):
    poly, px, py = generateRandomPolynomial(np.random.default_rng(seed), 7, 65536)
    y = [np.max(np.abs(poly(np.linspace(0, 1, k, endpoint=False)))) for k in nSample]

    if np.max(y) >= 1e6:
        print(seed)
        x = np.linspace(0, 1, 2**16)
        plt.plot(x, poly(x))
        plt.scatter(px, py, zorder=2)
        plt.grid()
        plt.show()

    result.append(y)

# result = np.diff(result, 1, 1)
# result = np.mean(result, 0)
result = np.max(result, 0)
plt.bar([str(n) for n in nSample], result, zorder=2)
# plt.bar([str(n) for n in nSample][:-1], result, zorder=2)
plt.grid()
plt.show()
