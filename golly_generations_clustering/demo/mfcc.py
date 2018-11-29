import argparse
import numpy
import python_speech_features
import sklearn.cluster
import sklearn.manifold
import shutil
import soundfile
from pathlib import Path

def extract_feature(path, n_frame=19):
    """
    mfcc.shape = (n_frame, n_cepstrum)
    """
    data, samplerate = soundfile.read(str(path))
    nfft = 1024
    mfcc = python_speech_features.mfcc(
        data,
        samplerate,
        winlen=nfft / samplerate,
        winstep=0.01,
        numcep=26,
        nfilt=52,
        nfft=nfft,
        preemph=0.97,
        ceplifter=22,
    )

    if mfcc.shape[0] < n_frame:
        zeros = numpy.zeros((n_frame - mfcc.shape[0], mfcc.shape[1]))
        mfcc = numpy.concatenate((mfcc, zeros), axis=0)
    elif mfcc.shape[0] > n_frame:
        mfcc = mfcc[0:n_frame]

    return numpy.ravel(mfcc)

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

    write_result(Path("cluster_mfcc"), n_clusters, cluster.labels_, filepath)

    numpy.save("data/filepath.npy", filepath)
    numpy.save("data/features.npy", features)
    numpy.save("data/features_shape.npy", (19, 26))
    numpy.save("data/labels.npy", cluster.labels_)
    numpy.save("data/centers.npy", cluster.cluster_centers_)
