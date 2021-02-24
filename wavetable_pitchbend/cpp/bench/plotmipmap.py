import matplotlib.pyplot as plt
import numpy as np
import soundfile
from pathlib import Path

import pprint
pp = pprint.PrettyPrinter(indent=2)

paths = {}
for path in Path("snd").glob("MipmapOscDeSoras*.wav"):
    name = path.stem
    name = name.replace("MipmapOscDeSoras", "")
    paths[int(name)] = path
paths = dict(sorted(paths.items(), reverse=True))
paths = {k: v for k, v in paths.items() if k >= 8}

fig, ax = plt.subplots(len(paths), 1)
fig.set_size_inches(6, 8)
fig.set_tight_layout(True)
for idx, path in enumerate(paths.values()):
    data, samplerate = soundfile.read(str(path))

    axis = ax[idx]
    lines = axis.stem(data, linefmt="black", markerfmt="r.")
    lines[0].zorder = 10
    for ln in lines[1:]:
        ln.set_linewidth(1)
        ln.set_color("black")
    axis.set_ylabel(f"{len(paths) - idx - 1}", rotation=0)
    axis.set_yticks([])
    axis.set_xticks([])
    axis.set_yticklabels([])
    axis.set_xticklabels([])
    axis.set_frame_on(False)
plt.show()
