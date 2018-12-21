import soundfile
from pathlib import Path

seconds = 0
for path in Path("snd/all").glob("*.wav"):
    print(f"\r\x1b[2K{str(path)}",end="")
    data, rate = soundfile.read(str(path))
    seconds += len(data) / rate
print(seconds)
