import numpy
import matplotlib.pyplot as pyplot

def bsplineN(i, k, x, t):
    if k == 1:
        if t[i] <= x and x < t[i + 1]:
            return 1
        else:
            return 0
    b1 = bsplineN(i, k - 1, x, t)
    if b1 != 0:
        b1 *= (x - t[i]) / (t[i + k - 1] - t[i])
    b2 = bsplineN(i + 1, k - 1, x, t)
    if b2 != 0:
        b2 *= (t[i + k] - x) / (t[i + k] - t[i + 1])
    return b1 + b2

def bsplineR(points, k, x, knot):
    return sum([p * bsplineN(i, k, x, knot) for i, p in enumerate(points)])

def plot_bspline():
    k = 4
    knot = numpy.hstack((
        numpy.zeros(k - 1),
        numpy.linspace(0, 1, k),
        numpy.ones(k - 1),
    ))
    points = [1, 2, 2, 1]
    x = numpy.linspace(0, 1, 1024)
    y = [bsplineR(points, k, D, knot) for D in x]
    pyplot.plot(x, y)

    px = numpy.linspace(0, 1, len(points))
    pyplot.scatter(px, points)

    pyplot.grid()
    pyplot.legend()
    pyplot.show()

def B0(t):
    return -t**3 / 6 + t**2 / 2 - t / 2 + 1 / 6

def B1(t):
    return t**3 / 2 - t**2 + 2 / 3

def B2(t):
    return -t**3 / 2 + t**2 / 2 + t / 2 + 1 / 6

def B3(t):
    return t**3 / 6

def polysum(t, p):
    return p[0] * B0(t) + p[1] * B1(t) + p[2] * B2(t) + p[3] * B3(t)

def J1B0(t):
    return -t**4 / 24 + t**3 / 6 - t**2 / 4 + t / 6 + J1B1(1)

def J1B1(t):
    return t**4 / 8 - t**3 / 3 + (2 * t) / 3 + J1B2(1)

def J1B2(t):
    return -t**4 / 8 + t**3 / 6 + t**2 / 4 + t / 6 + J1B3(1)

def J1B3(t):
    return t**4 / 24

def J2B0(t):
    return -t**5 / 120 + t**4 / 24 - t**3 / 12 + t**2 / 12 + (23 * t) / 24 + J2B1(1)

def J2B1(t):
    return t**5 / 40 - t**4 / 12 + t**2 / 3 + t / 2 + J2B2(1)

def J2B2(t):
    return -t**5 / 40 + t**4 / 24 + t**3 / 12 + t**2 / 12 + t / 24 + J2B3(1)

def J2B3(t):
    return t**5 / 120

def J2B0residue(t):
    return -t**5 / 120 + t**4 / 24 - t**3 / 12 + t**2 / 12 - t / 24 \
        + J2B1residue(1)

def J2B1residue(t):
    return t**5 / 40 - t**4 / 12 + t**2 / 3 - t / 2 + J2B2(1)

print(J2B1(1))

x = numpy.linspace(0, 1, 256)

# y = numpy.hstack((B3(x), B2(x), B1(x), B0(x)))
# pyplot.plot(y)

pyplot.title("Integrated Basis of Uniform Cubic B-Spline")
pyplot.xlabel("$t$")
# pyplot.ylabel("$B_i(t)$")

pyplot.plot(x + 0, B3(x), label="B3", ls="-", alpha=0.2, color="black")
pyplot.plot(x + 1, B2(x), label="B2", ls="--", alpha=0.4, color="black")
pyplot.plot(x + 2, B1(x), label="B1", ls="-", alpha=0.7, color="black")
pyplot.plot(x + 3, B0(x), label="B0", ls="--", alpha=1.0, color="black")

pyplot.plot(x + 0, J1B3(x), label="J B3", ls="-", alpha=0.2, color="blue")
pyplot.plot(x + 1, J1B2(x), label="J B2", ls="--", alpha=0.4, color="blue")
pyplot.plot(x + 2, J1B1(x), label="J B1", ls="-", alpha=0.7, color="blue")
pyplot.plot(x + 3, J1B0(x), label="J B0", ls="--", alpha=1.0, color="blue")

pyplot.plot(x + 0, J2B3(x), label="J^2 B3", ls="-", alpha=0.2, color="red")
pyplot.plot(x + 1, J2B2(x), label="J^2 B2", ls="--", alpha=0.4, color="red")
pyplot.plot(x + 2, J2B1(x), label="J^2 B1", ls="-", alpha=0.7, color="red")
pyplot.plot(x + 3, J2B0(x), label="J^2 B0", ls="--", alpha=1.0, color="red")

pyplot.plot(x + 0, J2B3(x), label="J^2 B3 residue", ls="-", alpha=0.2, color="orange")
pyplot.plot(x + 1, J2B2(x), label="J^2 B2 residue", ls="--", alpha=0.4, color="orange")
pyplot.plot(
    x + 2, J2B1residue(x), label="J^2 B1 residue", ls="-", alpha=0.7, color="orange")
pyplot.plot(
    x + 3, J2B0residue(x), label="J^2 B0 residue", ls="--", alpha=1.0, color="orange")

pyplot.legend(ncol=2)
pyplot.grid()
pyplot.show()
