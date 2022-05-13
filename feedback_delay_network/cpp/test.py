import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path

Path("img").mkdir(parents=True, exist_ok=True)

for path in Path("json").glob("*.json"):
    with open(path, "r", encoding="utf-8") as fi:
        data = json.load(fi)

    mat = np.array(data)
    dotted = mat.dot(mat.T)
    print(
        path.stem,
        mat.shape,
        f"det={np.linalg.det(mat)}",
    )
    # print(path.stem, mat, dotted, sep="\n\n", end="\n\n")

    cmap = plt.get_cmap("magma")
    plt.figure()
    plt.title(path.stem)
    for idx, row in enumerate(dotted):
        plt.plot(np.abs(row), color=cmap(idx / row.shape[0]), label=str(idx))
    plt.savefig(f"img/isunitary_{path.stem}.png")
    plt.close()

    plt.figure()
    plt.title(f"{path.stem}, Aboslute Values")
    CS = plt.imshow(np.abs(mat), cmap=cmap)
    plt.colorbar(CS)
    plt.savefig(f"img/absmatrix_{path.stem}.png")
    plt.close()

    plt.figure()
    plt.title(f"{path.stem}, Raw Values")
    CS = plt.imshow(mat, cmap=cmap)
    plt.colorbar(CS)
    plt.savefig(f"img/rawmatrix_{path.stem}.png")
    plt.close()
