import numpy
import matplotlib.pyplot as pyplot

def expPoly(α, β, time):
    return time**α * numpy.exp(-β * time)

def expPolyNormalized(α, β, time):
    normalize = expPoly(α, β, α / β)
    return time**α * numpy.exp(-β * time) / normalize

def plotAlpha(time):
    cmap = pyplot.get_cmap("viridis")
    pyplot.figure(figsize=(10, 5))
    alpha = numpy.geomspace(0.1, 16, 16)
    for idx, aa in enumerate(alpha):
        env = expPolyNormalized(aa, 1, time)
        pyplot.plot(
            time, env, alpha=0.9, lw=1, color=cmap(idx / len(alpha)), label=f"α={aa:.3f}")
    pyplot.title("Ê(t),  (β=1)")
    pyplot.xlabel("Time [s]")
    pyplot.ylabel("Amplitude")
    pyplot.grid()
    pyplot.legend(loc=1)
    pyplot.tight_layout()
    pyplot.show()

def plotBeta(time):
    cmap = pyplot.get_cmap("viridis")
    pyplot.figure(figsize=(10, 5))
    beta = numpy.geomspace(0.1, 16, 16)
    for idx, bb in enumerate(beta):
        env = expPolyNormalized(1, bb, time)
        pyplot.plot(
            time, env, alpha=0.9, lw=1, color=cmap(idx / len(beta)), label=f"β={bb:.3f}")
    pyplot.title("Ê(t),  (α=1)")
    pyplot.xlabel("Time [s]")
    pyplot.ylabel("Amplitude")
    pyplot.grid()
    pyplot.legend(loc=1)
    pyplot.tight_layout()
    pyplot.show()

samplerate = 48000
duration = 32

time = numpy.linspace(0, duration, duration * samplerate)

plotAlpha(time)
plotBeta(time)
