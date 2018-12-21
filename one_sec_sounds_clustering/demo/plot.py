import argparse
import functools
import itertools
import matplotlib.pyplot as pyplot
import multiprocessing
import numpy
import re
import soundfile
import subprocess
from pathlib import Path
from extract import Sound

class PlotData:
    def __init__(self, n_sound):
        self.index = 0
        self.n_frame = 100
        self.cmap_names = itertools.cycle(["plasma", "viridis", "cividis"])
        self.style_names = itertools.cycle(["dark_background"])
        self.samplerate = 44100
        self.duration_per_sound = 0.25

        self.wav_data = self.create_wav_data_buffer(n_sound)

    def get_start_position_in_wav_data(self):
        return int(self.index * self.samplerate * self.duration_per_sound)

    def create_wav_data_buffer(self, n_sound):
        length = int((n_sound * self.duration_per_sound + 1) * self.samplerate)
        return numpy.zeros(length)

    def write_wav(self, path):
        peak = numpy.max(numpy.abs(self.wav_data))
        soundfile.write(str(path), self.wav_data / peak, self.samplerate)

def reshape_mfcc_delta(data, x, y):
    if data is None:
        return (None, None)
    return numpy.hsplit(data.reshape(x, y).transpose(), 2)

def reshape_pitch(data, x, y):
    if data is None:
        return None
    return data.reshape(x, y).transpose()

def plot(plot_index, sound, center, label, cmap_name, n_frame):
    env_center, mfcc_center, pitch_center = center

    mfcc_numcep = 26
    n_pitch = 18

    mfcc_center, delta_center = reshape_mfcc_delta(mfcc_center, n_frame * 2,
                                                   mfcc_numcep)
    mfcc, delta = reshape_mfcc_delta(sound.mfccdelta, n_frame * 2, mfcc_numcep)

    pitch_center = reshape_pitch(pitch_center, n_frame, n_pitch)
    pitch = reshape_pitch(sound.pitch, n_frame, n_pitch)

    cmap = pyplot.get_cmap(cmap_name)

    fig = pyplot.figure(figsize=(19.20, 10.80), dpi=10)
    fig.suptitle(
        f"Cluster{label}, {sound.filepath.stem}",
        fontsize=16,
        x=0.05,
        horizontalalignment="left",
        color=cmap(0.5))
    axs = fig.subplots(2, 3)

    # top mid, mfcc center
    ax = axs[0][1]
    ax.set_title("MFCC of Cluster Center")
    ax.set_xlabel("Time Step [0.01s/step]")
    ax.set_ylabel("Cepstral Coefficent")
    ax.set_xlim([0, n_frame])
    ax.set_ylim([0, mfcc_numcep])
    if mfcc_center is None:
        ax.text(
            n_frame / 2,
            mfcc_numcep / 2,
            "N/A",
            fontsize=48,
            ha="center",
            va="center")
    else:
        ax.pcolormesh(mfcc_center, cmap=cmap)

    # bottom mid, mfcc sound
    ax = axs[1][1]
    ax.set_title("MFCC of Data Point")
    ax.set_xlabel("Time Step [0.01s/step]")
    ax.set_ylabel("Cepstral Coefficent")
    ax.set_xlim([0, n_frame])
    ax.set_ylim([0, mfcc_numcep])
    ax.pcolormesh(mfcc, cmap=cmap)

    # top right, delta center
    ax = axs[0][2]
    ax.set_title("Delta Cluster Center")
    ax.set_xlabel("Time Step [0.01s/step]")
    ax.set_ylabel("Cepstral Coefficent")
    ax.set_xlim([0, n_frame])
    ax.set_ylim([0, mfcc_numcep])
    if delta_center is None:
        ax.text(
            n_frame / 2,
            mfcc_numcep / 2,
            "N/A",
            fontsize=48,
            ha="center",
            va="center")
    else:
        ax.pcolormesh(delta_center, cmap=cmap)

    # bottom right, delta sound
    ax = axs[1][2]
    ax.set_title("Delta of Data Point")
    ax.set_xlabel("Time Step [0.01s/step]")
    ax.set_ylabel("Cepstral Coefficent")
    ax.set_xlim([0, n_frame])
    ax.set_ylim([0, mfcc_numcep])
    ax.pcolormesh(delta, cmap=cmap)

    # top left, envelope
    ax = axs[0][0]
    ax.set_title("Envelope")
    ax.set_xlabel("Time Step [0.01s/step]")
    ax.set_ylabel("Peak Value")
    ax.set_ylim([0, 1])
    if env_center is not None:
        ax.plot(env_center, label="cluster center", color=cmap(0.2))
    ax.plot(sound.envelope, label="data point", color=cmap(0.8))
    ax.legend()

    # bottom left, some pitch related feature
    ax = axs[1][0]
    ax.set_title("Pitch")
    ax.set_xlabel("Index")
    ax.set_ylabel("Cent")
    ax.set_ylim([-6000, 6000])
    if pitch_center is not None:
        ax.plot(pitch_center[0], label="cluster center", color=cmap(0.2))
        for pt in pitch_center[1:]:
            ax.plot(pt, color=cmap(0.2))
    ax.plot(pitch[0], label="data point", color=cmap(0.8))
    for pt in pitch[1:]:
        ax.plot(pt, color=cmap(0.8))
    ax.legend()

    fig.tight_layout()
    pyplot.subplots_adjust(top=0.92)

    fig.savefig(f"video/{plot_index:08d}.png", dpi=100)

    pyplot.close()

def open_dir(root, centers, sounds, plot_data, depth=0):
    center_path = root / Path("centers.npy")
    if depth < len(centers) and center_path.exists():
        centers[depth] = numpy.load(str(center_path))
        paths = sorted([path for path in root.iterdir() if path.is_dir()],
                       key=str)
        for path in paths:
            open_dir(path, centers, sounds, plot_data, depth + 1)
        return

    pattern = re.compile(r"(\d+$)|(outlier$)")
    center = [None for _ in range(len(centers))]
    indices = []
    for i in range(depth):
        try:
            index = pattern.search(root.parts[i - depth]).group(0)
        except AttributeError:
            print("Unexpected directory name at " + root)
            return

        indices.append(index)
        if index == "outlier":
            continue
        center[i] = (centers[i][int(index)])

    cmap_name = next(plot_data.cmap_names)
    pyplot.style.use(next(plot_data.style_names))
    label = "-".join(indices)

    print(f"\r\x1b[2K{indices}", end="")  # debug, "\x1b[2K" erase line.

    sound_plot = []
    start_index = plot_data.index
    for path in root.glob("*.wav"):
        sound = None
        for index, sound in enumerate(sounds):
            if sound.filepath.stem == path.stem:
                sound = sounds.pop(index)
                break
        sound_plot.append(sound)

        data, _ = soundfile.read(str(sound.filepath))
        start = plot_data.get_start_position_in_wav_data()
        end = start + len(data)
        plot_data.wav_data[start:end] += data
        plot_data.index += 1

    indices = [i for i in range(start_index, plot_data.index)]
    plot_partial = functools.partial(
        plot,
        center=center,
        label=label,
        cmap_name=cmap_name,
        n_frame=plot_data.n_frame,
    )
    plot_args = [(index, sound) for index, sound in zip(indices, sound_plot)]
    with multiprocessing.Pool() as pool:
        pool.starmap(plot_partial, plot_args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract features from wav file.")
    parser.add_argument(
        "src_name",
        metavar="src_name",
        type=str,
        help="source directory name. snd/src_name")
    args = parser.parse_args()

    Path("video").mkdir(parents=True, exist_ok=True)

    sounds = numpy.load(f"data/{args.src_name}.npy")

    plot_data = PlotData(len(sounds))

    root = Path("cluster") / Path(args.src_name)
    clusters = ["envelope", "mfcc", "pitch"]
    open_dir(
        root,
        [None for _ in range(len(clusters))],
        sounds.tolist(),
        plot_data,
    )
    print()

    plot_data.write_wav("video/merged.wav")

    subprocess.run([
        "ffmpeg",
        "-framerate",
        str(1 / plot_data.duration_per_sound),
        "-i",
        "video/%08d.png",
        "-i",
        "video/merged.wav",
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        "18",
        "-r",
        "30",
        "-c:a",
        "aac",
        "-b:a",
        "256k",
        "-y",
        "video/output.mp4",
    ])
