"""
ディレクトリ raw_recording 内のすべてのファイルを加工して、ディレクトリ snd/all に出力する。
"""

import numpy
import soundfile
import subprocess
from pathlib import Path

out_dir = Path("snd/all")
out_dir.mkdir(exist_ok=True)

temp_path = out_dir / Path("🚧temp🚧.wav")

for path in Path("raw_recording").iterdir():
    process = subprocess.run(["soxi", "-t", path], stdout=subprocess.DEVNULL)
    if process.returncode != 0:
        continue

    process = subprocess.run(
        ["sox", "--norm", path, temp_path, "remix", "-", "rate", "-v", "44.1k"])
    if process.returncode != 0:
        continue

    try:
        data, samplerate = soundfile.read(str(temp_path))
    except:
        print("Failed to read " + str(path))
        continue

    print(f"\r\x1b[2KSlicing {str(path)}", end="")

    digits = len(str(int(len(data) / samplerate)))
    index = 0
    start = 0
    end = samplerate  # 1秒間隔で切るので samplerate * 1。

    def write_wav(data):
        soundfile.write(
            str(out_dir / Path(f"{path.stem}_{index:0{digits}d}.wav")),
            data,
            samplerate,
        )

    while end < len(data):
        write_wav(data[start:end])
        index += 1
        start = end
        end += samplerate
    write_wav(data[start:])
print()

temp_path.unlink()
