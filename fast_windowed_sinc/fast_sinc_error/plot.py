import json
import matplotlib.pyplot as plt
import numpy as np


def plotForEach(data):
    cmap = plt.get_cmap("plasma")
    for key, value in data.items():
        plt.figure()
        plt.title(key)
        for idx, (cutoff, fir) in enumerate(value.items()):
            plt.plot(
                fir,
                color=cmap(idx / len(value)),
                alpha=0.75,
                lw=1,
                label=f"{float(cutoff):.3f}",
            )
        plt.legend()
        plt.grid()
        plt.tight_layout()
    plt.show()


def plotCompare(data):
    ref = data["ref"]
    fast = data["fast"]

    nCutoff = len(ref)
    nCol = int(np.ceil(np.sqrt(nCutoff)))
    nRow = nCol - (nCol * nCol - nCutoff) // nCol

    tap = len(ref[next(iter(ref))])
    xticks = [0, tap // 2, tap // 2 + 1, tap - 1]

    fig, ax = plt.subplots(nRow, nCol)
    fig.set_size_inches((16, 8))
    for idx, cutoff in enumerate(ref.keys()):
        axis = ax[idx // nCol][idx % nCol]
        axis.set_title(f"cutoff={float(cutoff):.3f}")
        axis.plot(ref[cutoff], color="blue", alpha=0.75, lw=1, label="ref.")
        axis.plot(fast[cutoff], color="red", alpha=0.75, lw=1, label="fast")
        axis.set_xticks(xticks, labels=[""] * len(xticks))
        if idx == 0:
            axis.legend()
        axis.grid()
    for idx in range(nCutoff, nCol * nRow):
        ax[idx // nCol][idx % nCol].axis("off")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    with open("test_windowed_sinc.json", "r", encoding="utf-8") as fp:
        data = json.load(fp)

    # plotForEach(data)
    plotCompare(data)
