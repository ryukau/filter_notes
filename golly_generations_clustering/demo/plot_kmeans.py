import itertools
import matplotlib
import matplotlib.pyplot as pyplot
import numpy

def get_errors(labels, features, centers):

    errors = numpy.zeros_like(centers)
    n_samples = numpy.zeros((centers.shape[0], 1), dtype=numpy.int32)

    for label, feature in zip(labels, features):
        errors[label] += numpy.abs(feature - centers[label])
        n_samples[label][0] += 1

    return (numpy.ravel(n_samples), errors / n_samples)

def plot(name, data, n_samples, winstep, x_len, y_len, row, col):
    data = data.reshape((data.shape[0], x_len, y_len)).transpose((0, 2, 1))

    fig = pyplot.figure(figsize=(12.80, 7.20), dpi=10)

    axs = fig.subplots(row, col)

    color_norm = matplotlib.colors.Normalize(
        numpy.min(data),
        numpy.max(data),
    )
    for c, r in itertools.product(range(col), range(row)):
        index = c + r * col
        if index >= len(data):
            break
        ax = axs[r][c]
        ax.set_title(f"{index}, n={n_samples[index]}")
        ax.axis("off")
        ax.pcolormesh(
            data[index],
            cmap=pyplot.get_cmap("magma"),
            norm=color_norm,
        )

    fig.tight_layout(pad=0, h_pad=0, w_pad=-2)
    fig.savefig(f"img/{name}.png", dpi=100)
    pyplot.close()

labels = numpy.load("data/labels.npy")
features = numpy.load("data/features.npy")
centers = numpy.load("data/centers.npy")
x_len, y_len = numpy.load("data/features_shape.npy")

n_samples, errors = get_errors(labels, features, centers)

plot("errors", errors, n_samples, 0.01, x_len, y_len, 5, 8)
plot("centers", centers, n_samples, 0.01, x_len, y_len, 5, 8)
