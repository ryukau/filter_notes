import argparse
import numpy
import scipy.signal
import sklearn.cluster
import shutil
import soundfile
from pathlib import Path

def extract_feature(path, n_frame=39):
    """
    スペクトログラムは2次元のデータ。
    sklearn.cluster で使えるように numpy.ravel で1次元にして返す。

    spectrogram.shape = (n_freq, n_frame)

    n_frame のデフォルト値はテストに使ったデータセットを調べて決めた。
    """
    data, samplerate = soundfile.read(str(path))
    frequency, time, spectrogram = scipy.signal.spectrogram(data, samplerate)

    if spectrogram.shape[1] < n_frame:
        zeros = numpy.zeros((spectrogram.shape[0],
                             n_frame - spectrogram.shape[1]))
        spectrogram = numpy.concatenate((spectrogram, zeros), axis=1)
    elif spectrogram.shape[1] > n_frame:
        spectrogram = spectrogram[:][0:n_frame]

    max_value = numpy.max(spectrogram)
    min_value = numpy.min(spectrogram)

    # return numpy.ravel(numpy.transpose(spectrogram))
    # return numpy.ravel(numpy.transpose((spectrogram - min_value) / max_value))

    spectrogram = numpy.log((spectrogram - min_value) / max_value)
    spectrogram[spectrogram == -numpy.inf] = -2000
    return numpy.ravel(numpy.transpose(spectrogram))

def write_result(output_directory, n_clusters, labels, filepath):
    if output_directory.exists():
        shutil.rmtree(output_directory)

    digits = len(str(abs(n_clusters - 1)))
    output_directories = [
        output_directory / Path(f"{index:0{digits}d}")
        for index in range(n_clusters)
    ]

    for directory in output_directories:
        directory.mkdir(parents=True, exist_ok=True)

    for label, path in zip(labels, filepath):
        shutil.copy(path, output_directories[label])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract features from wav file.")
    parser.add_argument(
        "src_path", metavar="src_path", type=str, help="Source directory path.")
    args = parser.parse_args()

    directory_path = Path(args.src_path)
    if not directory_path.is_dir():
        print("Invalid path.")

    filepath = [path for path in directory_path.glob("*.wav")]
    features = numpy.array([extract_feature(path) for path in filepath])

    n_clusters = 40
    cluster = sklearn.cluster.KMeans(n_clusters=n_clusters).fit(features)

    write_result(
        Path("cluster_spectrogram"), n_clusters, cluster.labels_, filepath)

    numpy.save("data/filepath.npy", filepath)
    numpy.save("data/features.npy", features)
    numpy.save("data/features_shape.npy", (39, 129))
    numpy.save("data/labels.npy", cluster.labels_)
    numpy.save("data/centers.npy", cluster.cluster_centers_)
