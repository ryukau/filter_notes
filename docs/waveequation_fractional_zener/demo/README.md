# Fractional Zener Model Wave Equation
このディレクトリには [fractional Zener model wave equation](http://heim.ifi.uio.no/~sverre/papers/2011_HolmNasholm-fractZener-JournAcoustSocAm.pdf) のシミュレーションを行う3つのデモがあります。

`draw.py` でGUI版のシミュレーションを実行できます。かなり重たいです。

```bash
python3 draw.py
```

`./render.sh` で連番のpng画像と mp4/h264 の動画をimgディレクトリに書き出します。動画の書き出しには `ffmpeg` が必要です。

```bash
./render.sh
```

# 使用ライブラリ
- PyQt5
- NumPy

# License
MIT
