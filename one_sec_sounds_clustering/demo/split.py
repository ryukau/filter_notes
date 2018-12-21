"""
ディレクトリ raw_synth 内のすべてのファイルを加工して、ディレクトリ snd/all に出力する。
"""

import subprocess
from pathlib import Path

out_dir = Path("snd/all")
out_dir.mkdir(parents=True, exist_ok=True)
for path in Path("raw_synth").glob("*"):
    channel = subprocess.check_output(
        ['soxi', '-c', str(path)],
        encoding="utf-8",
    ).rstrip()
    for ch in range(1, int(channel) + 1):
        out_name = out_dir / Path(f"{path.stem}_ch{ch}_normed.wav")
        proc = subprocess.run([
            "sox", "--norm", path, out_name, "remix",
            str(ch), "rate", "-v", "44.1k"
        ])
