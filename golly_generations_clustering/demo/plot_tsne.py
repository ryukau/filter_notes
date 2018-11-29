import matplotlib.pyplot as pyplot
import numpy

def sort_points_by_label(labels, embedding, n_clusters):
    bag = {n: [] for n in range(n_clusters)}
    for label, point in zip(labels, embedding):
        bag[label].append(point)
    return bag

labels = numpy.load("data/labels.npy")
embedding = numpy.load("data/embedding.npy")

n_clusters = numpy.max(labels) + 1

points_dict = sort_points_by_label(labels, embedding, n_clusters)

cmap = pyplot.get_cmap("nipy_spectral")
palette = [cmap(i) for i in numpy.linspace(0, 1, n_clusters)]

pyplot.figure(figsize=(19.20, 10.80), dpi=10)
for label, points in points_dict.items():
    axis = numpy.transpose(numpy.array(points))
    pyplot.scatter(
        axis[0],
        axis[1],
        s=48,
        c=palette[label],
        marker=f"${label}$",
        linewidths=None,
        edgecolors="none",
        label=str(label))
pyplot.savefig(f"img/tsne.png", dpi=100)
# pyplot.show()
