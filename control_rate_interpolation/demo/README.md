実行には次のライブラリが必要です。

- Python3
  - [NumPy](https://numpy.org/)
  - [matplotlib](https://matplotlib.org/)
  - [SoundFile](https://pysoundfile.readthedocs.io/en/latest/)
- C++17
  - [libsndfile](http://www.mega-nerd.com/libsndfile/)

## Windows
`run.ps1` を実行するとビルドします。 libsndfile を vcpkg からリンクします。

```ps1
./run.ps1
```

## Bash 環境
`run.sh` を実行するとビルドします。 `g++` が必要です。

```bash
./run.sh
```
