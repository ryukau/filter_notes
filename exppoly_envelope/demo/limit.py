import numpy
import matplotlib.pyplot as pyplot

def process(dtype):
    fmax = numpy.finfo(dtype).max
    alpha = dtype(128)

    timeMax = numpy.power(fmax, dtype(1) / alpha, dtype=dtype)
    print(timeMax)

    samplerate = dtype(48000)
    duration = dtype(32)
    time = numpy.linspace(dtype(0), duration, duration * samplerate, dtype=dtype)
    t_α = numpy.power(time, alpha, dtype=dtype)
    t_α[numpy.isinf(t_α)] = 0

    pyplot.plot(time, t_α / fmax)
    pyplot.axvline(timeMax, color="red")
    # pyplot.xlim([timeMax - 0.1, timeMax + 0.1])
    pyplot.show()

# process(numpy.float32)
# process(numpy.float64)

def printTimeMax(dtype):
    fmax = numpy.finfo(dtype).max
    alpha = numpy.array([2**i for i in range(1, 16)], dtype=dtype)
    timeMax = numpy.power(fmax, dtype(1) / alpha, dtype=dtype)
    print(dtype)
    for a, t in zip(alpha, timeMax):
        print(f"{a:9.1f}, {t}")

printTimeMax(numpy.float32)
printTimeMax(numpy.float64)
