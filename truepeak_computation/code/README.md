# トゥルーピークの計算の調査に使用したコード
以下のリンクから本記事を読めます。

- [トゥルーピークの計算を読む (github.io)](../truepeak_computation.html)

## 使用ライブラリ
実効には以下のライブラリが必要です。

### C++
- [libsndfile](http://www.mega-nerd.com/libsndfile/)
- [nlohmann/json: JSON for Modern C++](https://github.com/nlohmann/json)

### Python 3
- [SciPy](https://www.scipy.org/)
- [NumPy](https://numpy.org/)
- [Matplotlib](https://matplotlib.org/)
- [PySoundFile](https://pysoundfile.readthedocs.io/en/latest/)
- [CVXOPT](http://cvxopt.org/)

## コードの実行
### SOCP FIR の設計
`fractionaldelaysocp.py` を実行すると [SOCP FIR](https://ryukau.github.io/filter_notes/fractional_delay_filter_socp/fractional_delay_filter_socp.html) が設計できます。また `plotresponse.py` を実行すると設計した SOCP FIR の周波数、位相、群遅延特性がプロットできます。

```bash
python fractionaldelaysocp.py
python plotresponse.py
```

実行後に出力される `socp.json` に FIR フィルタの係数が入っています。プロットには `socp.json` が `plotresponse.py` と同じディレクトリに配置されている必要がありおます。またコマンドラインに C++ の `std::array` で使える形式にフォーマットしたフィルタ係数が出力されます。

SOCP FIR のパラメータを変えるときは `fractionaldelaysocp.py` を開いて以下の行を変更します。

```python
if __name__ == "__main__":
    table = createTable(7, 5, omega_max=0.65, omega_density=1) # ここでパラメータを設定。
```

### Thiran オールパスフィルタの設計
`thiran_allpass` ディレクトリ以下に関連するコードをまとめています。

設計するには `thiran.py` を実行します。

```bash
python thiran.py
```

コマンドラインに C++ の `std::array` で使えるようにフォーマットされた、 [sos](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.sosfilt.html) 形式のフィルタ係数が出力されます。

Thiran オールパスフィルタのパラメータを設定するには `thiran.py` を開いて `printSos` の引数を変更してください。

```python
printSos(order=8, oversample=4) # パラメータ設定例。
```

`thiran.wxmx` は [wxMaxima](https://wxmaxima-developers.github.io/wxmaxima/) で開くことができる Maxima のコードです。

### EBU TECH 3341 のテスト
以下のリンクから EBU TECH 3341 のテスト信号がダウンロードできます。

- [EBU TECH 3341, 3342, 3343 のテスト信号のダウンロードページ](https://tech.ebu.ch/publications/ebu_loudness_test_set)

ダウンロードした zip を解凍して `code/data` に配置します。初回は `data` ディレクトリを作成してください。 以下の `tree` コマンドの出力のように、 `data/ebu-loudness-test-setv05` 内に wav ファイルが入っていれば `ebutest.py` が実行できます。

```
+---data
|   |
|   +---ebu-loudness-test-setv05
|   |       1kHz Sine -20 LUFS-16bit.wav
|   |       1kHz Sine -26 LUFS-16bit.wav
|   |       1kHz Sine -40 LUFS-16bit.wav
|   |       EBU-reference_listening_signal_pinknoise_500Hz_2kHz_R128.wav

# ...
```

### 計算速度の比較
計算速度の比較では `data` ディレクトリ以下に配置された音声ファイルを再帰的にすべて読み込みます。読めないファイルは自動的にスキップします。

libsndfile で読めるファイル形式にだけ対応しています。 [libsndfile のトップページ](http://www.mega-nerd.com/libsndfile/)に対応ファイル形式の表があります。

以下のリンク先に freesound.org から入手した、実験に使ったデータの一覧を掲載しています。

- [データセット一覧を見る](data/dataset.md)

ダウンロードを終えたら `data` ディレクトリ内に解凍してください。

あとは C++ のコードをコンパイルして実行すればベンチマークが始まります。以下は Windows 10 で vcpkg を使っている場合のコマンドです。

```bash
cd cpp/benchmark
mkdir build
cd build

cmake -DCMAKE_TOOLCHAIN_FILE="/src/vcpkg/scripts/buildsystems/vcpkg.cmake" ..
cmake --build . --config Release
./Release/bench.exe
```

**注意:** `bench.exe` は `build` ディレクトリ内から実行しないと `data` ディレクトリを読み取りません。実験用のコードなので、手間を省くために汎用性のある実装にしていません。気になるときは `measure` 関数の `dirPath` あたりを変えてみてください。

Linux の場合はを以下のコマンドに変更してみてください。

```bash
cd cpp/benchmark
mkdir build
cd build

cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
./Release/bench
```

**注意:** Linux 環境で動作確認していないので `cpp/benchmark/CMakeLists.txt` の変更が必要になるかもしれません。

### トゥルーピークの近似値の誤差の測定
計算速度の比較と同様に `data` ディレクトリ以下に適当な音声データを配置します。後は以下の順番で Python 3 のスクリプトを実行してください。

```
python measuresinc.py # Sinc 補間によるトゥルーピークの真値を計算。
python measuresocp.py # SOCP FIR によるトゥルーピークの近似値を計算。
python error.py       # 誤差の測定とプロット。
```

**注意:** 32 bit の Python 3 では `measuresinc.py` の実行時にメモリ不足になることがあります。 64 bit の Python 3 を使ってください。

**注意:** `measuresinc.py` と `measuresocp.py` は実行終了までに数十分から数時間ほどの長い時間がかかります。

以下のファイルが出力されます。

- `measure_sinc.json`: Sinc 補間によるトゥルーピークの真値。
- `measure_socp.json`: SOCP FIR によるトゥルーピークの近似値。
- `measure_socp_filter.json`: 測定した SOCP FIR のフィルタ係数。

### BS.1770 FIR より効率のいいパラメータ
計算速度の比較と同様に `data` ディレクトリ以下に適当な音声データを配置します。後は以下の順番で Python 3 のスクリプトを実行してください。

```
python measuresinglesocp.py
python socpebutest.py
```

以下のファイルが出力されます。

- `error_socp_omega_max.json`: オーバーサンプリングの倍率を固定して `omega_max` を変えたときの SOCP FIR の誤差。
- `socp_ebutest_status.json`: EBU TECH 3341 のトゥルーピークのテスト通過の成否。 `true` ならテスト通過。

### 最悪の場合の信号
`worstsignal.py` を実行するとトゥルーピークが最大になる、最悪の場合の信号をレンダリングします。

```bash
python worstsignal.py
```

`worstsinc.py` は有限の長さの信号のトゥルーピークの最大値を計算できます。そのまま実行すると 1 秒、 1 分、 1 時間、 1 日の長さの信号のトゥルーピークの最大値を表示します。

```bash
python worstsinc.py
```

`cpp/sincerror` に `worstsinc.py` と同じ値を for 文で愚直に計算する実装があります。 Windows 10 で `cl.exe` がインストールされた環境では以下のコマンドで実行できます。 `cl.exe` は Visual Studio の C++ 環境をインストールするとパスが通っているはずです。

```
cd cpp/sincerror
cl /O2 /EHsc /std:c++17 .\sincerror.cpp
.\sincerror.exe
```

Linux では以下を試してみてください (動作未確認) 。

```
cd cpp/sincerror
g++ -O3 -std=c++17 sincerror.cpp # clang++ でも可。
./a.out
```

## ライセンス
MIT
