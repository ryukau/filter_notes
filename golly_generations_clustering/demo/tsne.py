import numpy
import sklearn.manifold

features = numpy.load("data/features.npy")
tsne = sklearn.manifold.TSNE(n_components=2).fit(features)
numpy.save("data/embedding.npy", tsne.embedding_)
