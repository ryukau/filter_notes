import argparse
import collections
import functools
import itertools
import matplotlib.pyplot as pyplot
import numpy
import sklearn.cluster
import shutil
from pathlib import Path
from extract import Sound

class Order:
    def __init__(self, func, feature, params):
        self.func = func
        self.feature = feature
        self.params = params

def kmeans(features, n_clusters):
    cluster = sklearn.cluster.KMeans(
        n_clusters=n_clusters,
        # verbose=1,
    ).fit(features)
    return (n_clusters, cluster.labels_, cluster.cluster_centers_)

def affinity_propagation(features, damping):
    cluster = sklearn.cluster.AffinityPropagation(
        damping=damping,
        max_iter=2000,
        convergence_iter=30,
        # verbose=True,
    ).fit(features)
    return (len(cluster.cluster_centers_), cluster.labels_,
            cluster.cluster_centers_)

def compose_func(func, points, params):
    if params is not None:
        return functools.partial(func, points, *params)
    return functools.partial(func, points)

def sort_data_by_label(data, n_clusters, labels):
    bag = [[] for _ in range(n_clusters)]
    for dat, label in zip(data, labels):
        bag[label].append(dat)
    return bag

def recurse_clustering(depth, orders, sounds):
    if len(sounds) <= 1:
        return sounds

    order = orders[depth]

    clustering_func = compose_func(
        order.func,
        [getattr(sound, order.feature) for sound in sounds],
        order.params,
    )

    n_clusters, labels, centers = clustering_func()
    sounds = sort_data_by_label(sounds, n_clusters, labels)

    depth += 1
    if depth >= len(orders):
        sounds.append(centers)
        return sounds

    sounds = [recurse_clustering(depth, orders, sound) for sound in sounds]
    sounds.append(centers)
    return sounds

def start_clustering(orders, sounds):
    return recurse_clustering(0, orders, sounds)

def delete_previous_result(data_prefix):
    output_directory = Path("cluster") / Path(data_prefix)
    if output_directory.exists():
        shutil.rmtree(output_directory)

def get_next_path(n_data, name, index, digits):
    if n_data > 1:
        return Path(f"{name}{index:0{digits}d}")
    return Path(f"{name}_outlier")

def write_result(data_prefix, sounds, feature_names, depth=0):
    """
    sounds = [
        [[s0, s1, ...], centers],
        [[s0, s1, ...], centers],
        ...
        [[s0, s1, ...], centers],
        centers,
    ]
    """
    output_directory = Path("cluster") / data_prefix
    output_directory.mkdir(parents=True, exist_ok=True)

    if isinstance(sounds[-1], numpy.ndarray):
        numpy.save(str(output_directory / "centers.npy"), sounds.pop())

    if isinstance(sounds[0], list):
        digits = len(str(abs(len(sounds) - 1)))
        for index, sound in enumerate(sounds):
            next_prefix = data_prefix / get_next_path(
                len(sound), feature_names[depth], index, digits)
            write_result(next_prefix, sound, feature_names, depth + 1)
        return

    for path in [snd.filepath for snd in sounds]:
        shutil.copy(path, output_directory)

def print_clusters(clusters, depth=0):
    space = " " * depth
    for c in clusters:
        if isinstance(c, list):
            print(space + "[",)
            print_clusters(c, depth + 1)
            print(space + "],")
        elif isinstance(c, Sound):
            pass
        else:
            print(space + str(type(c)) + ",")

def clustering(data_prefix):
    sounds = numpy.load(f"data/{data_prefix}.npy")

    orders = [
        Order(affinity_propagation, "envelope", (0.9,)),
        Order(affinity_propagation, "mfccdelta", (0.6,)),
        Order(affinity_propagation, "pitch", (0.5,)),
    ]

    delete_previous_result(data_prefix)

    # print_clusters(start_clustering(orders, sounds))
    write_result(
        Path(data_prefix),
        start_clustering(orders, sounds),
        [order.feature for order in orders],
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract features from wav file.")
    parser.add_argument(
        "src_name",
        metavar="src_name",
        type=str,
        help="source directory name. snd/src_name")
    args = parser.parse_args()

    clustering(args.src_name)
