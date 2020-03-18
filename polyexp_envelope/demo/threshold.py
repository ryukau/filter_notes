import numpy
import matplotlib.pyplot as pyplot
import scipy.special.lambertw as lambertw

def getTime(α, β, x, normalize, k=0):
    return -α * lambertw(-β * (x * normalize)**(1 / α) / α, k) / β

def polyExp(α, β, time):
    return time**α * numpy.exp(-β * time)

def polyExpNormalized(α, β, time):
    normalize = polyExp(α, β, α / β)
    return (
        time**α * numpy.exp(-β * time) / normalize,
        normalize,
    )

samplerate = 48000
duration = 8

time = numpy.linspace(0, duration, duration * samplerate)

alpha = 1
beta = 1

curve, normalize = polyExpNormalized(alpha, beta, time)

xx = numpy.linspace(1, 0, 10)
# xx = 1e-5
tt0 = getTime(alpha, beta, xx, normalize, 0)
value_tt0, _ = polyExpNormalized(alpha, beta, tt0)
tt1 = getTime(alpha, beta, xx, normalize, -1)
value_tt1, _ = polyExpNormalized(alpha, beta, tt1)
tt2 = getTime(alpha, beta, xx, normalize, 1)
value_tt2, _ = polyExpNormalized(alpha, beta, tt2)

pyplot.plot(time, curve)
pyplot.scatter(tt0, value_tt0, color="red", label="k=0")
pyplot.scatter(tt1, value_tt1, color="orange", label="k=-1")
pyplot.scatter(tt2, value_tt2, color="purple", label="k=1")
pyplot.title("Get time from amplitude of Ê(t)")
pyplot.xlabel("Time [s]")
pyplot.ylabel("Amplitude")
pyplot.grid()
pyplot.legend()
pyplot.show()
