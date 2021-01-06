# ピークホールドによるエンベロープ
音のリミッタで使うために、ピークホールドを使ったエンベロープの検出を実装します。

## 前向きホールド
入力信号の絶対値を計算します。

<figure>
<img src="img/input_signal.svg" alt="Plot of input signal and its absolute value." style="padding-bottom: 12px;"/>
</figure>

絶対値を前から順に読み取って、そのときまでの最大値を出力します。これがもっとも単純な前向きホールドです。

<figure>
<img src="img/simple_forward_peak_hold.svg" alt="Plot of simple forward peak hold." style="padding-bottom: 12px;"/>
</figure>

しかし、単純な前向きホールドでは減衰するエンベロープをうまく検出できません。

<figure>
<img src="img/simple_forward_peak_hold_of_decaying_signal.svg" alt="Plot of simple forward peak hold of decaying signal." style="padding-bottom: 12px;"/>
</figure>

そこで、一定時間が経ったらホールドしていた最大値を忘れるようにします。この文章で前向きホールドと言うときは、最大値を忘れる実装を指すことにします。

<figure>
<img src="img/forward_peak_hold_of_decaying_signal.svg" alt="Plot of forward peak hold." style="padding-bottom: 12px;"/>
</figure>

Python 3 で実装します。この文章の Python 3 のコードは上から順にテキストファイルにコピペすると動きます。実行には [SciPy](https://www.scipy.org/) 、 [NumPy](https://numpy.org/) 、 [matplotlib](https://matplotlib.org/) が必要です。

```python
import matplotlib.pyplot as plt
import numpy as np

def peakHoldForward(sig, holdtime, reset=0):
    """
    holdtime の単位はサンプル数。
    最後に最大値が更新されてから holdtime サンプル後に hold の値をリセットする。
    """
    out = np.empty(len(sig))
    hold = reset
    counter = 0
    for i in range(len(sig)):
        if counter > 0:  # 注意: 比較に >= を使うとホールド時間が 1 サンプル延びる。
            counter -= 1
        else:
            hold = reset
        if hold <= sig[i]:  # 注意: カウンタをリセットしたいので <= で比較。
            hold = sig[i]
            counter = holdtime
        out[i] = hold
    return out

rng = np.random.default_rng(0)

sig = rng.normal(0, 1, 256)
absed = np.abs(sig)
holdF = peakHoldForward(absed, 32, 0)

plt.plot(sig, label="Input")
plt.plot(absed, label="Absed")
plt.plot(holdF, label="Forward Hold")
plt.show()
```

## スムーシング
ある信号の絶対値をなだらかに包み込むような信号のことを、その信号のエンベロープと言います。しきい値やニーなどの特性を加えたエンベロープの逆数を入力信号に掛け合わせることで音量を制限するリミッタが作れます。

前向きホールドによるエンベロープをそのまま使うとホールド値を忘れるときにポップノイズが出てしまうので、フィルタを使ってスムーシングします。ここで使うフィルタは、一定時間同じ値が入力されると、その値に到達することが保証されているものとします。移動平均フィルタの重ね掛けや、 Bessel フィルタはこのような設計ができます。リミッタ向きのフィルタについては[「ステップ応答が S 字を描くフィルタ」](../s_curve_step_response_filter/s_curve_step_response_filter.html)にまとめています。

以降では、以下のコードで設計した三角窓の FIR フィルタを使います。

```python
def makeTriangleFir(delay):
    """delay は 2 より大きい整数。"""
    fir = np.interp(
        np.linspace(0, 1, delay + 1),
        [0, 0.5, 1],
        [0, 1, 0],
    )
    return fir / np.sum(fir)
```

前向きホールドの出力にフィルタをかけると立ち上がりと立ち下りが鈍るので、入力信号との間でピークの位置がずれます。以降ではホールドの出力にフィルタをかけた信号をエンベロープと呼びます。

<figure>
<img src="img/smoothed_hold.svg" alt="Plot of smoothed hold." style="padding-bottom: 12px;"/>
</figure>

そこで入力信号にディレイをかけてエンベロープとピークの位置をあわせます。実装によりますが、リミッタのレイテンシはこのディレイから来ています。

<figure>
<img src="img/smoothed_hold_and_delayed_signal_arrow.svg" alt="Plot of smoothed hold and delayed signal." style="padding-bottom: 12px;"/>
</figure>

ここで致命的な問題があります。上の図の赤い矢印で示した箇所などで、入力信号のピークがエンベロープからはみ出しています。このようなはみ出しが起こるエンベロープをリミッタで使うと音量が制限できません。

以下は、はみ出しが起こるパターンを示した図です。入力信号のピークはフィルタをかける前のエンベロープからは、はみ出しません。しかし、ディレイを加えた入力信号はフィルタをかけたエンベロープからはみ出すことがあります。

<figure>
<img src="img/protruding_peak_fail.svg" alt="Image of a pattern which causes protruding peak." style="padding-bottom: 12px;"/>
</figure>

はみ出しが起こるのは、ホールド中に入力された小さいピークを見逃しているからです。下の図のように、ホールド値を忘れるときに入力を振り返って見逃したピークをホールドしなおせば、はみ出しが防げそうです。

<figure>
<img src="img/protruding_peak_desired.svg" alt="Image of desired envelope which eliminates protruding peak." style="padding-bottom: 12px;"/>
</figure>

## 理想的なホールド
ホールド値より小さいピークを見逃さないためには、ホールド時間と同じ長さのバッファを用意して入力を覚えておく方法があります。出力にはバッファ内の最大値を使います。以下はこのアルゴリズムのコード例です。

```python
from collections import deque

def idealPeakHoldNaive(sig, holdtime: int):
    out = np.empty(len(sig))
    buffer = deque([0 for _ in range(holdtime)])
    counter = 0
    for i in range(len(sig)):
        buffer.append(abs(sig[i]))
        buffer.popleft()
        out[i] = max(buffer)
    return out
```

ここでは上のコードのアルゴリズムのことを理想的なホールドの素朴な実装と呼びます。素朴な実装の問題点は計算量です。上のコードの `out[i] = max(buffer)` は Python では 1 行で書けますが、計算量はバッファの長さに比例します。したがってホールド時間を長くすると計算が重たくなります。

### リアルタイムで使えるアルゴリズム {#realtime-implementation}
計算量の問題を解決するために試行錯誤したところ、局所最大値をキューに保存することで高速に計算するアルゴリズムを思いつきました。

例として以下のような入力について考えます。図の赤い実線が理想的なホールドの出力です。赤い点線は出力値の候補を表しています。

<figure>
<img src="img/ideal_hold_local_maximum.svg" alt="Plot of local maximum of ideal hold. 1st peak has the largest value. 2nd peak is the smallest. The value of 3rd peak is in between 1st and 2nd." style="padding-bottom: 12px;"/>
</figure>

ここでピーク 2 は以下の条件が揃っているので出力に影響を及ぼしません。

- ピーク 2 はピーク 1 のホールド時間内に現れる。
- ピーク 2 より後に、より値の大きいピーク 3 が現れる。
- ピーク 3 もピーク 1 のホールド時間内に現れる。
- ピーク 2 はピーク 1 とピーク 3 よりも小さい。

言い換えると、ピーク 3 が入力されたときにホールド時間だけ入力を振り返って、ピーク 3 よりも小さい値を忘れてしまってもいいということです。

#### レシピ
以下のメモリを確保します。

- ホールド時間 + 1 サンプルのディレイ。 0 で初期化。
- ホールド時間のサンプル数と同じ長さのキュー。 0 で初期化。

計算手順です。

1. キューの最後尾が入力サンプル `x0` 以上の値になるまで、キューの最後尾を繰り返し除去。
2. `x0` をキューの後ろに追加。
3. `x0` をディレイに入力。
4. もしディレイの出力とキューの先頭の値が等しければ、キューから値を取り出す。
5. キューに一つでも値が残っていれば、キューの先頭を出力。そうでなければ 0 を出力。

上のレシピに沿った Python 3 による実装です。

```python
from collections import deque

def idealPeakHoldFast(sig, holdTime):
    """
    sig: 入力信号。
    holdTime: ホールド時間。単位はサンプル数。
    """
    out = np.empty(len(sig))                     # 出力信号。
    buffer = deque([0 for _ in range(holdTime)]) # ディレイのバッファ。
    hold = deque([])                             # ホールド値のキュー。
    for i in range(len(sig)):
        x0 = sig[i]

        if len(hold) > 0: # 手順 1
            idx = len(hold) - 1
            while idx >= 0:
                if hold[idx] < x0:
                    hold.pop()
                else:
                    break
                idx -= 1

        hold.append(x0)   # 手順 2

        buffer.append(x0) # 手順 3
        delayOut = buffer.popleft()

        if len(hold) > 0 and delayOut == hold[0]: # 手順 4
            hold.popleft()

        out[i] = hold[0] if len(hold) > 0 else 0 # 手順 5

    return out
```

以下は理想的なホールドの出力です。高速な実装と理想的な実装の出力が一致しています。

<figure>
<img src="img/ideal_hold.svg" alt="Plot of ideal hold envelope." style="padding-bottom: 12px;"/>
</figure>

スムーシングした理想的なホールドの出力です。ディレイをかけた入力がエンベロープからはみ出していないことが確認できます。

<figure>
<img src="img/smoothed_ideal_hold.svg" alt="Plot of smoothed hold envelope." style="padding-bottom: 12px;"/>
</figure>

### テスト
理想的なホールドを使えばはみ出しが起こらないことをテストします。

テスト用の入力として、ランダムにインパルスが散らばったノイズを作ります。

```python
def nextTime(rng, rate):
    """
    ポアソン過程。
    もし単位時間当たり N 回のイベントが起こるなら `rate` を `1/N` に設定する。
    """
    return -np.log(1.0 - rng.uniform(0, 1)) / rate

def pulseNoise(rng, rate, length):
    out = np.zeros(length)
    time = 0
    while time < length:
        time += nextTime(rng, rate)

        t1 = int(time)
        t2 = t1 + 1
        amp = rng.uniform(0, 1)
        if t1 < length:
            out[t1] = amp * (t2 - time)
        if t2 < length:
            out[t2] = amp * (time - t1)
    return out
```

ポアソン過程の実装は Jeff Preshing さんによる以下の記事を参考にしました。

- [How to Generate Random Timings for a Poisson Process](https://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/)

以下はテストコードです。

```python
from multiprocessing import Pool

def verifyHoldEnvelope(sig, env, name):
    condition = (sig - env) > 8 * np.finfo(np.float64).eps
    if np.any(condition):
        idx = np.where(condition)
        diff = np.take(sig, idx) - np.take(env, idx)
        print(f"Test failed: {name}", diff, sep="\n")

def jobTestIdealHoldRandom(seed):
    hold = 32   # ホールド時間 (サンプル数)
    delay = 32  # フィルタの遅延 (サンプル数)
    triangleFir = makeTriangleFir(delay)

    rng = np.random.default_rng(seed)
    sig = pulseNoise(rng, 1 / 8, 48000)

    naive = idealPeakHoldNaive(sig, hold)
    fast = idealPeakHoldFast(sig, hold)

    smoothedNaive = signal.lfilter(triangleFir, 1, naive)
    smoothedFast = signal.lfilter(triangleFir, 1, fast)

    delayed = np.hstack((np.zeros(hold), sig[:-hold]))

    verifyHoldEnvelope(delayed, smoothedNaive, f"Naive, {seed}, {delay}")
    verifyHoldEnvelope(delayed, smoothedFast, f"Fast, {seed}, {delay}")

def testIdealHoldRandom(nTest=1024):
    with Pool() as pool:
        seeds = [i for i in range(nTest)]
        idx = 0
        for result in pool.imap_unordered(jobTestIdealHoldRandom, seeds):
            print(f"\r{idx}", end="")  # 処理が終わったシード値を表示。
            idx += 1

if __name__ == "__main__":
    testIdealHoldRandom()
```

フィルタの遅延はホールド時間以下の任意のサンプル数に設定できますが、値によっては誤差が増えます。上のコードでは誤差の上限を適当に `8 * np.finfo(np.float64).eps` としています。

はみ出しを検知したときは `Test failed: ...` を出力します。シード値 0 から 1023 まで今回実装した理想的なホールドをテストしたところ、はみ出しは検知されませんでした。

### C++ による実装
以下のコードの `PeakHold` が理想的なホールドです。

```c++
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

// 整数サンプルのディレイ。
template<typename Sample> class Delay {
public:
  std::vector<Sample> buf;
  size_t wptr = 0;
  size_t rptr = 0;

  Delay(size_t size = 65536) : buf(size) {}

  void resize(size_t size)
  {
    buf.resize(size);
    wptr = 0;
    rptr = 0;
  }

  void reset() { std::fill(buf.begin(), buf.end(), Sample(0)); }

  void setFrames(size_t delayFrames)
  {
    if (delayFrames >= buf.size()) delayFrames = buf.size();
    rptr = wptr - delayFrames;
    if (rptr >= buf.size()) rptr += buf.size(); // Unsigned overflow case.
  }

  Sample process(Sample input)
  {
    if (++wptr >= buf.size()) wptr -= buf.size();
    buf[wptr] = input;

    if (++rptr >= buf.size()) rptr -= buf.size();
    return buf[rptr];
  }
};

// メモリの確保と解放を減らした std::deque の代用データ構造。
template<typename T> struct RingQueue {
  std::vector<T> buf;
  size_t wptr = 0;
  size_t rptr = 0;

  void resize(size_t size) { buf.resize(size); }

  void reset(T value = 0) { std::fill(buf.begin(), buf.end(), value); }

  inline size_t size()
  {
    auto sz = wptr - rptr;
    if (sz >= buf.size()) sz += buf.size(); // Unsigned overflow case.
    return sz;
  }

  inline bool empty() { return wptr == rptr; }

  T &front() { return buf[increment(rptr)]; }
  T &back() { return buf[wptr]; }

  inline size_t increment(size_t idx)
  {
    if (++idx >= buf.size()) idx -= buf.size();
    return idx;
  }

  inline size_t decrement(size_t idx)
  {
    if (--idx >= buf.size()) idx += buf.size(); // Unsigned overflow case.
    return idx;
  }

  void push_back(T value)
  {
    wptr = increment(wptr);
    buf[wptr] = value;
  }

  T pop_front()
  {
    rptr = increment(rptr);
    return buf[rptr];
  }

  T pop_back()
  {
    wptr = decrement(wptr);
    return buf[wptr];
  }
};

/*
理想的なホールドの実装。
- `setFrames(0)` のとき出力はすべて 0 。
- `setFrames(1)` のとき入力をバイパス。
*/
template<typename Sample> struct PeakHold {
  Sample neutral = 0;
  Delay<Sample> delay;
  RingQueue<Sample> hold;

  PeakHold(size_t size = 65536) {
    resize(size);
    setFrames(1);
  }

  void resize(size_t size)
  {
    delay.resize(size);
    hold.resize(size);
  }

  void reset()
  {
    delay.reset();
    hold.reset(neutral);
  }

  void setFrames(size_t frames) { delay.setFrames(frames); }

  Sample process(Sample x0)
  {
    if (!hold.empty()) {
      for (size_t idx = hold.size(); idx > 0; --idx) {
        if (hold.back() < x0)
          hold.pop_back();
        else
          break;
      }
    }

    hold.push_back(x0);

    auto delayOut = delay.process(x0);
    if (!hold.empty() && delayOut == hold.front()) hold.pop_front();

    return !hold.empty() ? hold.front() : neutral;
  }
};

int main()
{
  // PeakHold の使用例。
  constexpr size_t length = 64;

  PeakHold<float> hold;
  hold.setFrames(4);

  std::minstd_rand rng(0);
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  for (size_t idx = 0; idx < length; ++idx) {
    auto input = dist(rng);
    std::cout << input << ", " << hold.process(input) << "\n";
  }

  return 0;
}
```

`PeakHold` ではニュートラル値 `neutral` を変更できるようにしています。

`RingQueue` は `push_back` や `pop_front` でメモリの確保や解放が行われないようにした `std::deque` の代替です。速度はほとんど同じです。

## その他
### 後ろ向きホールド
試行錯誤しているときに下の図のような後ろ向きホールドを思いついたのですが、理想的なホールドの生成には使えないことがわかりました。

<figure>
<img src="img/backward_hold.svg" alt="Example plot of backward hold." style="padding-bottom: 12px;"/>
</figure>

前向きホールドと組み合わせて最大値を取ることで理想的なホールドに近づけられます。入力信号を $x$ 、ホールド時間を $d$ 、理想的なホールドを $H_I$ 、 前向きホールドを $H_F$ 、後ろ向きホールドを $H_B$ とすると以下の計算です。

$$
H_I(x) = \max(H_F(x), H_B(x) z^{-d})
$$

$z^{-d}$ は $d$ サンプルの遅延を表しています。

後ろ向きホールドの実装です。インデックスの並びが逆転している以外は、前向きホールドと同じです。

```python
def peakHoldBackward(sig, holdtime, reset=0):
    out = np.empty(len(sig))
    hold = reset
    counter = 0
    for i in reversed(range(len(sig))): # ここだけ変更。
        if counter > 0:
            counter -= 1
        else:
            hold = reset
        if hold <= sig[i]:
            hold = sig[i]
            counter = holdtime
        out[i] = hold
    return out
```

以下は前向きと後ろ向きホールドを組み合わせても理想的なホールドの生成に失敗するケースです。黒い縦線が入力、青い縦線がディレイをかけた入力です。

<figure>
<img src="img/backward_hold_fail.svg" alt="Example plots of a case that backward hold fails to make correct envelope." style="padding-bottom: 12px;"/>
</figure>

前から順に a, b, c の3つピークがあります。以下の条件が揃うと b を見逃してしまうので失敗します。

- b は a のホールド時間内に現れる。
- c は b のホールド時間内に現れる。
- a > b かつ c > b 。
