# DR-110風のシンバル
コードを実行するには[Python3](https://www.python.org/)と以下のライブラリが必要です。

- [SciPy](https://www.scipy.org/)
- [NumPy](http://www.numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [pyFFTW](https://hgomersall.github.io/pyFFTW/)
- [PySoundFile](http://pysoundfile.readthedocs.io/en/0.9.0/)

以下のコマンドで `dr110_cymbal.py` を実行するとディレクトリ `snd` にDR-110風のシンバルがレンダリングされます。

```sh
python3 dr110_cymbal.py
```

`bandpass.py` と `metal_filter.py` はフィルタの特性をプロットするプログラムです。
