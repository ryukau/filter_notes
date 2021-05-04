import matplotlib.pyplot as plt
import numpy as np
import soundfile
from pathlib import Path

img_dir = Path("img")
if not img_dir.exists():
    img_dir.mkdir(parents=True, exist_ok=True)

for path in Path("snd").glob("*.wav"):
    data, samplerate = soundfile.read(str(path))

    plt.figure(figsize=(10, 8), tight_layout=True)
    plt.title(path.stem)
    plt.plot(np.arange(len(data)) / samplerate, data, lw=1, color="black")
    plt.grid()
    plt.savefig(str(img_dir / f"{path.stem}.svg"))
plt.show()
