# プロットギャラリー
クラスタリングの特徴として試した、素朴なスペクトログラム、正規化したスペクトログラム、MFCCから得られたプロットです。

素朴なスペクトログラムは `scipy.signal.spectrogram` をそのまま特徴に使っています。

正規化したスペクトログラムはデータごとに `scipy.signal.spectrogram` の範囲を [0, 1] に正規化した特徴です。

MFCCは `python_speech_features.mfcc` を使っています。

プロットは耳以外での評価が何かできないかと思ったので作りました。どれも曖昧でこれといった正解はありません。

## Elbow Method
見た目はどれも似たように見えますが、素朴なスペクトログラムだけエラーの値がやたら低いです。

### 素朴なスペクトログラム
<figure>
<img src="img/spectrogram_naive/elbow_method.png" alt="Image of elbow method plot from naive spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### 正規化したスペクトログラム
<figure>
<img src="img/spectrogram_normalized/elbow_method.png" alt="Image of elbow method plot from normalized spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### 対数スペクトログラム
<figure>
<img src="img/spectrogram_log/elbow_method.png" alt="Image of elbow method plot from log spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### MFCC
<figure>
<img src="img/mfcc/elbow_method.png" alt="Image of elbow method plot from mfcc features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

## t-SNE
スペクトログラムから得られたクラスタはどれもいまいちな感じがします。特徴の次元が高いときは直接K-Meansに渡すよりもt-SNEで次元を減らしてからK-Meansという手法も考えられます。

対数スペクトログラムは他と雰囲気が異なっています。

### 素朴なスペクトログラム
<figure>
<img src="img/spectrogram_naive/tsne.png" alt="Image of t-SNE plot from naive spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### 正規化したスペクトログラム
<figure>
<img src="img/spectrogram_normalized/tsne.png" alt="Image of t-SNE plot from normalized spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### 対数スペクトログラム
<figure>
<img src="img/spectrogram_log/tsne.png" alt="Image of t-SNE plot from log spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### MFCC
<figure>
<img src="img/mfcc/tsne.png" alt="Image of t-SNE plot from mfcc features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

## 中央値
中央値が表すスペクトログラムとMFCCです。各画像の縦が周波数で下から上に向かって大きくなります。横は時間で左から右に向かって進んでいます。明るい部分ほど係数が大きくなります。 `n` はクラスタに含まれるデータポイントの数を表しています。

スペクトログラムによるクラスタはどれも `n=1` となる outlier がMFCCに比べると多いです。3つの内では正規化したスペクトログラムが一番良さそうに見えます。対数スペクトログラムは極端な値が出ているように見えます。素朴なスペクトログラムが真っ黒なのは一部の値が極端に大きいからです。

### 素朴なスペクトログラム
<figure>
<img src="img/spectrogram_naive/centers.png" alt="Image of plot of K-Means centers from naive spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### 正規化したスペクトログラム
<figure>
<img src="img/spectrogram_normalized/centers.png" alt="Image of plot of K-Means centers from normalized spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### 対数スペクトログラム
<figure>
<img src="img/spectrogram_log/centers.png" alt="Image of plot of K-Means centers from log spectrogram features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### MFCC
<figure>
<img src="img/mfcc/centers.png" alt="Image of plot of K-Means centers from mfcc features." style="width: 600px;padding-bottom: 12px;"/>
</figure>

## 平均絶対誤差 (Mean Absolute Error)
中央値と、中央値が代表するクラスタ内のデータポイントとの平均絶対誤差です。

### 素朴なスペクトログラム

<figure>
<img src="img/spectrogram_naive/errors.png" alt="Image of plot of mean absolute errors between K-Means centers and corresponding data points from naive spectrogram." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### 正規化したスペクトログラム

<figure>
<img src="img/spectrogram_normalized/errors.png" alt="Image of plot of mean absolute errors between K-Means centers and corresponding data points from normalized spectrogram." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### 対数スペクトログラム

<figure>
<img src="img/spectrogram_log/errors.png" alt="Image of plot of mean absolute errors between K-Means centers and corresponding data points from log spectrogram." style="width: 600px;padding-bottom: 12px;"/>
</figure>

### MFCC

<figure>
<img src="img/mfcc/errors.png" alt="Image of plot of mean absolute errors between K-Means centers and corresponding data points from MFCC." style="width: 600px;padding-bottom: 12px;"/>
</figure>

## 耳で聞いた印象
### 素朴なスペクトログラム
`n<=20` のクラスタは似たような音が入っている。 `n` が大きいクラスタは何でもあり。

### 正規化したスペクトログラム

- 0 はいろいろ混ざっていた。
- 1 は大半がプチノイズでたまにカラーノイズ。
- 16 はいろいろ混ざっていた。0よりもトーンが高め。
- 3, 4, 5, 7, 13, 18, 20, 22, 36, 38 はカラーノイズ。
- 25, 27, 28, 30 は低めのトーン。
- 2, 10, 11, 14, 24 は高めのトーン。

### 対数スペクトログラム
- 0 はプチノイズのみ。
- 1 は何でもあり。
- 16 はプチノイズと高いトーン。
- 25 はカラーノイズとトーンが混じった音。ノイズのほうが強め。
- 39 はトーンが強めの音。

### MFCC
- 0, 2, 15, 20, 25, 28 は低めのトーン。
- 1, 11, 21 はプチノイズ。
- 4, 7, 16 はカラーノイズ。
- 3, 8, 9, 12, 13, 14, 19, 29, 31, 35, 39 は高いトーン。
- 17, 18, 30, 34, 36 は高めのカラーノイズと高めのトーンの中間のような音。
- 22, 23, 24, 26, 37, 38 はプチノイズと高いトーン。
- 32, 33 は低めのカラーノイズと低めのトーンの中間のような音。
