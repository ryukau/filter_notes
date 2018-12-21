# 1秒以下の音のクラスタリング
[GollyのGenerationsルールから生成した音のクラスタリング](../golly_generations_clustering/golly_generations_clustering.html)が意外とうまく行ったので、1秒以下の音のクラスタリングを試しました。

## コードについて
Python3では次のライブラリを使っています。

- [matplotlib](https://matplotlib.org/)
- [PySoundFile](https://pysoundfile.readthedocs.io/en/0.9.0/)
- [python_speech_features](https://python-speech-features.readthedocs.io/en/latest/#)
- [scikit-learn](https://scikit-learn.org/stable/index.html)
- [SciPy, NumPy](https://www.scipy.org/)

## データセット
データセットに含まれる音は、長さを1秒以下、サンプリング周波数を44100Hzにそろえています。

今回作ったデータセットの大きさは6016サンプルで長さは約81分です。

### シンセサイザの音
以前作ったシンセサイザから、それぞれ1000ほどのサンプルをレンダリングしました。

今回使ったシンセサイザです。

- [Singen0.2](https://ryukau.github.io/Singen0.2/)
- [Singen0.3](https://ryukau.github.io/Singen0.3/)
- [Pluck](https://ryukau.github.io/Pluck/)
- [KSCymbal](https://ryukau.github.io/KSCymbal/)
- [PADcymbal](https://ryukau.github.io/PADcymbal/)
- [FDNCymbal](https://ryukau.github.io/FDNCymbal/)
- [WaveCymbal](https://ryukau.github.io/WaveCymbal/)

ブラウザのデベロッパツールのコンソールから次のコードを実行してレンダリングしました。

```javascript
var title = document.getElementsByTagName("title")
var button = title[0].innerText === "Singen0.3"
  ? ui.buttonRandom : buttonRandom
var counter = 0
var id = window.setInterval(
  () => {
    if (counter > 1024) clearInterval(id)
    button.onClick()
    ++counter
  },
  1000
)
```

Pluck, KSCymbal, PADcymbal, FDNCymbal, WaveCymbalはステレオでレンダリングした音をチャンネルごとに分離してそれぞれ別のサンプルとして扱っています。

得られた音を次のコードで正規化しました。[SoX](http://sox.sourceforge.net/)を使っています。

- [シンセサイザの音を加工するコードを読む (github.com)](demo/split.py)

### フィールドレコーディングの音
近所を歩いてフィールドレコーディングした音をデータセットに加えました。

フィールドレコーディングの音はSoXの `remix -` でモノラルにしてから1秒間隔で切っています。

- [フィールドレコーディングの音を加工するコードを読む (github.com)](demo/slice.py)

- [SoX man page](http://sox.sourceforge.net/sox.html) - EFFECTSの節に `remix` の解説

## 特徴抽出
### MFCC+デルタ
`python_speech_features.mfcc` で取り出したMFCCを `numpy.ravel` で1次元にして `sklearn.cluster.AffinityPropagation` でクラスタリングしました。

今回使った `python_speech_features.mfcc` のパラメータです。

```python
import python_speech_features
import soundfile

data, samplerate = soundfile.read("path/to/wav_file")

nfft = 1024
mfcc = python_speech_features.mfcc(
    data,
    samplerate,
    winlen=nfft / samplerate,
    winstep=0.01,
    numcep=26,
    nfilt=52,
    nfft=nfft,
    lowfreq=0,
    highfreq=20000,
    preemph=0.0,
    ceplifter=0,
)
```

今回は1秒以下の音のクラスタリングなので MFCC のフレーム数が `1 / winstep = 100` になるように調整しました。フレーム数が100より少ないときは0で埋めたフレームを付け足しています。フレーム数が100より大きいときは、音の始まりから100フレームだけを取り出して、残りのフレームを切り捨てています。

[Practical Cryptography の MFCC チュートリアル](http://www.practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/)にMFCCのデルタで結果が改善することがあると書いてあったので `numpy.concatenate` でMFCCと、MFCCのデルタをつないで一つのデータポイントにまとめました。

```python
# mfcc の取得は省略。
delta = python_speech_features.base.delta(mfcc, 1)
mfccdelta = numpy.concatenate((numpy.ravel(mfcc), numpy.ravel(delta)))
```

ここでは `mfccdelta` のことをMFCC+デルタと呼んでいます。

### エンベロープ
ここでのエンベロープは信号を短い区間に区切って、それぞれの区間で信号の絶対値の最大値を取り出したものです。

```python
def get_envelope(data, samplerate, winstep, n_frame):
    """
    :data: 一次元の信号。
    :samplerate: サンプリング周波数。
    :winstep: 区間の長さ。秒。
    :n_frame: 区間の数。
    """
    data_abs = numpy.abs(data)
    envelope = numpy.zeros(n_frame)
    index = None
    step = int(winstep * samplerate)
    start = 0
    end = step
    for frame in range(n_frame):
        if end >= len(data_abs):
            index = frame
            break
        envelope[frame] = numpy.max(data_abs[start:end])
        start = end
        end += step
    if index is not None:
        envelope[index] = numpy.max(data_abs[start:])
        index += 1
    return envelope
```

エンベロープの区間の長さと数はMFCCの対応するパラメータと合わせました。区間の長さは `winstep=0.01` 、区間の数は `n_frame=100` です。

### ソートしたピッチ
ピッチは `python_speech_features.sigproc.framesig` で区切った各フレームから次の手順で取り出しました。

1. NSD type II の局所最大点を全て取り出す。
2. 局所最大点を降順にソート。
3. ソートした局所最大点の前から18個のインデックスを取り出して周波数を計算。
4. 周波数をセント値に変換。
5. セント値を昇順にソート。

手順 5. のソートがないと似たような音でまとまりにくくなります。推定したピッチは似たような値で順番が入れ替わっていることが多いので、ソートなしだとクラスタリングで計算される距離が大きくなることが予想されます。

CMND type II も試したのですが、テストに使った小さなデータセットでは NSD type II のほうが良い結果が出ました。

係数 k の値ではテストデータでの結果は変わりませんでした。

## クラスタリング
エンベロープ、MFCC+デルタ、ソートしたピッチは別物なのでクラスタリングを分けることにしました。以下はエンベロープ -> MFCC+デルタ -> ソートしたピッチの順でクラスタリングした結果の一例です。クラスタリングの階層構造がなんとなく見て取れるかと思います。 Outlier と判断されたクラスタは番号が飛んでいます。

```
cluster/test2
├── envelope1
│   ├── mfccdelta2
│   │   ├── pitch1
│   │   │   ├── fast_vib2.wav
│   │   │   ├── high_vib.wav
│   │   │   ├── klang4.wav
│   │   │   ├── noisy1.wav
│   │   │   └── vib2.wav
│   │   └── pitch_outlier
│   │       └── fast_vib1.wav
│   └── mfccdelta_outlier
│       ├── ping1.wav
│       └── ping2.wav
├── envelope2
│   ├── mfccdelta0
│   │   ├── pitch1
│   │   │   ├── ding.wav
│   │   │   ├── fm.wav
│   │   │   └── noisy4.wav
│   │   └── pitch_outlier
│   │       └── sweep_to_high.wav
│   └── mfccdelta_outlier
│       └── vib1.wav
├── envelope3
│   ├── mfccdelta0
│   │   ├── pitch1
│   │   │   ├── klang1.wav
│   │   │   ├── klang3.wav
│   │   │   └── noisy3.wav
│   │   └── pitch_outlier
│   │       ├── mid_crack.wav
│   │       └── noisy5.wav
│   └── mfccdelta_outlier
│       ├── klang2.wav
│       └── vib3.wav
└── envelope_outlier
    └── noisy2.wav
```

クラスタリング手法は `sklearn.cluster.AffinityPropagation` を使いました。エンベロープのクラスタリングでは `damping=0.9` 、MFCC+デルタのクラスタリングでは `damping=0.6` 、ソートしたピッチのクラスタリングでは `damping=0.5` としました。

次の動画はクラスタリングの結果です。音は0.25秒間隔で再生されます。残りの0.75秒は次の音と重なっています。

<iframe width="640" height="360" src="https://www.youtube.com/embed/4SleDuWeYT4" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

多少は似たような音が集まっている気がします。

## その他
ラベルのついていないデータセットでのクラスタリングは無謀です。

今回使ったエンベロープとMFCC+デルタはそういう無茶でもそれなりになんとかしてくれました。エンベロープとMFCC+デルタによるクラスタリングではビブラートのような細かいピッチの揺れやカラーノイズはうまく区別できていないように感じたのでピッチ推定に注目したのですがあまり大きく改善したようには感じませんでした。

## The Infinite Drum Machine の手法
Google AI Experiments の [The Infinite Drum Machine](https://experiments.withgoogle.com/drum-machine) で使われていた手法を参考にしたので紹介します。 The Infinite Drum Machine は [Kyle Mcdonald さんの AudioNotebooks](https://github.com/kylemcdonald/AudioNotebooks)
 の視覚化です。クラスタリングはインターフェイス上で音を表す点の色付けに使っているだけのようです。

1. 音の始まりから250msを切り取って [STFT](https://en.wikipedia.org/wiki/Short-time_Fourier_transform)。
2. STFTの結果を [0, 1] の範囲に正規化。
3. [`skimage.measure.block_reduce`](http://scikit-image.org/docs/dev/api/skimage.measure.html#skimage.measure.block_reduce) で正規化したSTFTの結果を縮小。
4. 縮小したSTFTの結果をデータポイントとして t-SNE で2次元にマッピング。
5. マッピングした空間で K-Means を使ってクラスタリング。

`block_reduce` は画像の特徴抽出で使われる関数です。STFTの結果は2次元なので、画像とみなして処理できるというのは面白いです。

# 参考サイト

- [Practical Cryptography](http://www.practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/)
- [GitHub - kylemcdonald/AudioNotebooks: Collection of notebooks and scripts related to audio processing and machine learning.](https://github.com/kylemcdonald/AudioNotebooks)
- [The Infinite Drum Machine by Manny Tan & Kyle McDonald | Experiments with Google](https://experiments.withgoogle.com/drum-machine)
- [GitHub - googlecreativelab/aiexperiments-drum-machine: Thousands of everyday sounds, organized using machine learning.](https://github.com/googlecreativelab/aiexperiments-drum-machine)
- [An Intuitive Explanation of Convolutional Neural Networks – the data science blog](https://ujjwalkarn.me/2016/08/11/intuitive-explanation-convnets/)
