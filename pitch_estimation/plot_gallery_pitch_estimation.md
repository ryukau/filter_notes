# ピッチ推定のプロットギャラリー

## サイン波
サイン波の周波数ごとの平均絶対誤差です。この形式のグラフの縦軸は対数スケールにしています。

<figure>
<img src="img/error_per_frequency_sin_wave.png" alt="Image of plot of mean absolute error for each sine wave signal." style="width: 640px;padding-bottom: 12px;"/>
</figure>

誤差が最大となった YIN-CMND type I は低周波域で大きな誤差が出ています。

6000Hzより上の誤差を除外した場合が気になったのでプロットしました。

<figure>
<img src="img/error_sin_wave_under_6000hz.png" alt="Image of plot of mean absolute error to sine wave signal under 6000Hz." style="width: 640px;padding-bottom: 12px;"/>
</figure>

MPM-CMND type II の誤差が最小となりました。サンプリング周波数の3/4よりも高い周波数ではピッチ推定が失敗しやすいのかもしれません。

## ノイズを加えたサイン波
サイン波の周波数ごとの平均絶対誤差です。

<figure>
<img src="img/error_per_frequency_sin_with_noise.png" alt="Image of plot of mean absolute error for each frequency of sin with noise signal." style="width: 640px;padding-bottom: 12px;"/>
</figure>

ノイズの比率ごとの平均絶対誤差です。

<figure>
<img src="img/error_per_noise_ratio_sin_with_noise.png" alt="Image of plot of mean absolute error for each noise_ratio of sin with noise signal." style="width: 640px;padding-bottom: 12px;"/>
</figure>

### 手法ごとの誤差
明るいほど誤差が大きく、白い部分は `nan` を表しています。グラフの右にあるカラーバーの数字を $a$ とすると誤差の値は $10^a$ です。

NSD は信号の大きさ1に対してノイズの大きさが3を超えるとうまく計算できないようです。

<figure>
<img src="img/error_sin_with_noise_yin_cmnd_type1.png" alt="Image of plot of mean absolute error to sin with noise signal. Method is yin_cmnd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_with_noise_yin_cmnd_type2.png" alt="Image of plot of mean absolute error to sin with noise signal. Method is yin_cmnd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_with_noise_yin_nsd_type1.png" alt="Image of plot of mean absolute error to sin with noise signal. Method is yin_nsd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_with_noise_yin_nsd_type2.png" alt="Image of plot of mean absolute error to sin with noise signal. Method is yin_nsd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_with_noise_mpm_cmnd_type1.png" alt="Image of plot of mean absolute error to sin with noise signal. Method is mpm_cmnd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_with_noise_mpm_cmnd_type2.png" alt="Image of plot of mean absolute error to sin with noise signal. Method is mpm_cmnd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_with_noise_mpm_nsd_type1.png" alt="Image of plot of mean absolute error to sin with noise signal. Method is mpm_nsd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_with_noise_mpm_nsd_type2.png" alt="Image of plot of mean absolute error to sin with noise signal. Method is mpm_nsd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

## AM変調
キャリアの周波数ごとの平均絶対誤差です。

<figure>
<img src="img/error_per_carrier_frequency_sin_am.png" alt="Image of plot of mean absolute error for each carrier frequency of sin AM signal." style="width: 640px;padding-bottom: 12px;"/>
</figure>

モジュレータの周波数ごとの平均絶対誤差です。

<figure>
<img src="img/error_per_modulator_frequency_sin_am.png" alt="Image of plot of mean absolute error for each modulator frequency of sin AM signal." style="width: 640px;padding-bottom: 12px;"/>
</figure>

### 手法ごとの誤差
明るいほど誤差が大きく、白い部分は `nan` を表しています。グラフの右にあるカラーバーの数字を $a$ とすると誤差の値は $10^a$ です。

<figure>
<img src="img/error_sin_am_yin_cmnd_type1.png" alt="Image of plot of mean absolute error to sin AM signal. Method is yin_cmnd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_am_yin_cmnd_type2.png" alt="Image of plot of mean absolute error to sin AM signal. Method is yin_cmnd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_am_yin_nsd_type1.png" alt="Image of plot of mean absolute error to sin AM signal. Method is yin_nsd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_am_yin_nsd_type2.png" alt="Image of plot of mean absolute error to sin AM signal. Method is yin_nsd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_am_mpm_cmnd_type1.png" alt="Image of plot of mean absolute error to sin AM signal. Method is mpm_cmnd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_am_mpm_cmnd_type2.png" alt="Image of plot of mean absolute error to sin AM signal. Method is mpm_cmnd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_am_mpm_nsd_type1.png" alt="Image of plot of mean absolute error to sin AM signal. Method is mpm_nsd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_am_mpm_nsd_type2.png" alt="Image of plot of mean absolute error to sin AM signal. Method is mpm_nsd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

## FM変調
キャリアの周波数ごとの平均絶対誤差です。

<figure>
<img src="img/error_per_carrier_frequency_sin_fm.png" alt="Image of plot of mean absolute error for each carrier frequency of sin FM signal." style="width: 640px;padding-bottom: 12px;"/>
</figure>

モジュレータの周波数ごとの平均絶対誤差です。

<figure>
<img src="img/error_per_modulator_frequency_sin_fm.png" alt="Image of plot of mean absolute error for each modulator frequency of sin FM signal." style="width: 640px;padding-bottom: 12px;"/>
</figure>

### 手法ごとの誤差
明るいほど誤差が大きく、白い部分は `nan` を表しています。グラフの右にあるカラーバーの数字を $a$ とすると誤差の値は $10^a$ です。

<figure>
<img src="img/error_sin_fm_yin_cmnd_type1.png" alt="Image of plot of mean absolute error to sin FM signal. Method is yin_cmnd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_fm_yin_cmnd_type2.png" alt="Image of plot of mean absolute error to sin FM signal. Method is yin_cmnd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_fm_yin_nsd_type1.png" alt="Image of plot of mean absolute error to sin FM signal. Method is yin_nsd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_fm_yin_nsd_type2.png" alt="Image of plot of mean absolute error to sin FM signal. Method is yin_nsd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_fm_mpm_cmnd_type1.png" alt="Image of plot of mean absolute error to sin FM signal. Method is mpm_cmnd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_fm_mpm_cmnd_type2.png" alt="Image of plot of mean absolute error to sin FM signal. Method is mpm_cmnd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_fm_mpm_nsd_type1.png" alt="Image of plot of mean absolute error to sin FM signal. Method is mpm_nsd_type1." style="width: 640px;padding-bottom: 12px;"/>
</figure>

<figure>
<img src="img/error_sin_fm_mpm_nsd_type2.png" alt="Image of plot of mean absolute error to sin FM signal. Method is mpm_nsd_type2." style="width: 640px;padding-bottom: 12px;"/>
</figure>
