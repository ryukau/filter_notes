import matplotlib.pyplot as pyplot
import numpy
import python_speech_features
import sklearn.cluster
import soundfile
from pathlib import Path

def get_inertia(n_clusters, features):
    cluster = sklearn.cluster.KMeans(n_clusters=n_clusters).fit(features)
    return cluster.inertia_

if __name__ == "__main__":
    features = numpy.load("data/features.npy")

    x_range = (2, 101)
    errors = [get_inertia(k, features) for k in range(*x_range)]
    k_value = [k for k in range(*x_range)]

    pyplot.figure(figsize=(12.80, 7.20), dpi=10)
    pyplot.plot(k_value, errors, lw=1, color="gray")
    pyplot.plot(k_value, errors, "o", markersize=3, color="black")
    pyplot.xlabel("n_clusters")
    pyplot.ylabel("error")
    pyplot.xlim((x_range[0], x_range[1] - 1))
    pyplot.grid()
    pyplot.savefig(f"img/elbow_method.png", dpi=100)
