<style>
body {
  max-width: 704px;
  margin: auto;
  padding: 32px 8px;
}

canvas {
  /* image-rendering: pixelated; */
  display: inline-block;
  border-style: solid;
  border-width: 0.2px;
  border-color: #000000;
}

.numberInputLabel {
  display: inline-block;
  text-align: right;
  width: 160px;
  margin: 0 8px 4px 0;
  padding-right: 8px;
  border-right: solid 3px #606060;
}

.numberInputRange {
  display: inline-block;
  max-width: 220px;
  width: 50%;
  margin: 0 8px 0 0;
}

.numberInputNumber {
  display: inline-block;
  max-width: 100px;
  width: 20%;
  margin: 0;
}

code {
  overflow-x: scroll;
  overflow-y: hidden;
  white-space: pre;
}

.katex {
  font-size: 1.3em !important;
}
</style>

# 分数階微積分の数値計算
[Loverroの "Fractional Calculus: History, Definitions and Applications for the
Engineer" ](https://web.archive.org/web/20051029113800/http://www.nd.edu/~msen/Teaching/UnderRes/FracCalc.pdf)に沿って、数値計算がしやすそうなGrünwald–Letnikovの[分数階微積分](https://en.wikipedia.org/wiki/Fractional_calculus)について見ていきます。

## 分数階微積分とは
ここではざっくりとした導入だけを行います。基本的な理論については次の資料が参考になります。

- [Fractional Calculus: History, Definitions and Applications for the Engineer](https://web.archive.org/web/20051029113800/http://www.nd.edu/~msen/Teaching/UnderRes/FracCalc.pdf)
- [Fractional Calculus - Xuru's Website](http://www.xuru.org/fc/TOC.asp)
- [Fractional order derivative of a function & fractional numbers' factorial. - YouTube](https://www.youtube.com/watch?v=7PCQmxlX9mU)
- [Fractional calculus - Wikipedia](https://en.wikipedia.org/wiki/Fractional_calculus)

分数階微分について考えます。微分演算子 $D$ を導入します。

$$
D f(x) = \frac{d}{dx} f(x)
$$

例えば $1 / 2$ 階微分を $2$ 回繰り返すと1階微分になります。

$$
D^{\frac{1}{2}} \left( D^{\frac{1}{2}} f(x) \right)
= D f(x)
$$

$1 / 3$ 階微分を $3$ 回繰り返すと1階微分になります。

$$
D^{\frac{1}{3}} \left(
D^{\frac{1}{3}} \left(
  D^{\frac{1}{3}} f(x) \right) \right)
= \frac{d}{dx} f(x)
$$

つまり $1 / q$ 階微分を $q$ 回繰り返すと1階微分になります。

$$
\underbrace{
  D^{\frac{1}{q}} \left( \ldots \left( D^{\frac{1}{q}} f(x) \right) \right)
}_\text{1/q differentiation for q times}

= \frac{d}{dx} f(x)
$$

分数階積分も同様に捉えることができます。積分演算子 $J$ を定義します。

$$
J f(x) = \int_{0}^{x} f(s) ds
$$

$1 / q$ 階積分を $q$ 回繰り返すと1階積分になります。

$$
\underbrace{
  J^{\frac{1}{q}} \left( \ldots \left( J^{\frac{1}{q}} f(x) \right) \right)
}_\text{1/q integration for q times}

= \int_{0}^{x} f(s) ds
$$

分数階微積分が物理的にどういう意味を持つのかはよくわかりませんが、[波動方程式への応用](http://heim.ifi.uio.no/~sverre/papers/2011_HolmNasholm-fractZener-JournAcoustSocAm.pdf)が提案されています。

## Grünwald–Letnikovの分数階微分
[$n$ 階微分の差分形](https://en.wikipedia.org/wiki/Differential_of_a_function#Higher-order_differentials)です。

$$
d^n f(x) = \lim_{h \to 0} \frac{1}{h^n}
\sum_{m=0}^{n} (-1)^{m} \binom{n}{m} f(x - mh)
\,,\qquad n \in \Bbb{N}
$$

$m > n$ のとき $\binom{n}{m} = 0$ なので総和の上限を無限にできます。

$$
d^n f(x) = \lim_{h \to 0} \frac{1}{h^n}
\sum_{m=0}^{\infty} (-1)^{m} \binom{n}{m} f(x - mh)
\,,\qquad n \in \Bbb{N}
$$

自然数 $n$ を実数 $\alpha$ に置き換えます。

$$
d^\alpha f(x) = \lim_{h \to 0} \frac{1}{h^\alpha}
\sum_{m=0}^{\infty} (-1)^{m} \binom{\alpha}{m} f(x - mh)
\,,\qquad \alpha \in \Bbb{R}
$$

この形がGrünwald–Letnikovの分数階微分です。

## 二項係数の計算
実数 $\alpha$ が入力される[二項係数](https://en.wikipedia.org/wiki/Binomial_coefficient) $\binom{\alpha}{m}$ を計算できる形にします。まずは自然数についての二項係数を展開します。

$$
\binom{n}{m}
= \frac{n!}{m!(n - m)!}
\,,\qquad n \in \Bbb{N}\,,\quad m \in \Bbb{N}
$$

階乗は[ガンマ関数](https://en.wikipedia.org/wiki/Gamma_function) $\Gamma(n)$ に置き換えられます。

$$
\Gamma(n) = (n - 1)!, \qquad n \in \Bbb{N}
$$

ガンマ関数を用いて二項係数を変形します。

$$
\binom{\alpha}{m}
= \frac{\Gamma(\alpha + 1)}{\Gamma(m + 1)\Gamma(\alpha - m + 1)}
\,,\qquad \alpha \in \Bbb{R}\,,\quad m \in \Bbb{R}
$$

Grünwald–Letnikovの分数階微分に代入します。

$$
d^{\alpha} f(x)
= \lim_{h \to 0} \frac{1}{h^{\alpha}} \sum_{m=0}^{\infty} (-1)^{m}
\frac{\Gamma(\alpha + 1)}{\Gamma(m + 1)\Gamma(\alpha - m + 1)} f(x - mh)
\,, \qquad \alpha \in \Bbb{R}
$$

これで実数 $\alpha$ が入力されても計算できる形になりました。

ガンマ関数に入力される値が $0$ あるいは負の整数のときは無限大になるので、計算できない場合があります。Grünwald–Letnikovの分数階微分の式では $\alpha > 0$ かつ $m > 0$ かつ $\alpha > m - 1$ のときは問題なく計算できます。

## ガンマ関数の計算
ガンマ関数の数値計算には[Lanczos approximation](https://en.wikipedia.org/wiki/Lanczos_approximation)が使えます。次のコードはJavaScriptでの実装です。

```javascript
function gamma(z) {
  var p = [
    676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012,
    9.9843695780195716e-6, 1.5056327351493116e-7
  ]

  if (z < 0.5) {
    return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z))
  }

  z -= 1
  var a = 0.99999999999980993
  for (var i = 0; i < p.length; ++i) {
    a += p[i] / (z + i + 1)
  }
  var t = z + p.length - 0.5
  return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * a
}
```

## Grünwald–Letnikovの分数階積分
$\alpha$ 階の分数階積分も計算できます。 $-\alpha$ 階微分の形ではガンマ関数に負の値が入ることになるので計算できない場合が出てきます。

$$
\begin{aligned}
d^{-\alpha} f(x)

&= \lim_{h \to 0} \frac{1}{h^{-\alpha}}
\sum_{m=0}^{\infty} (-1)^{m} \binom{-\alpha}{m} f(x - mh)\\

&= \lim_{h \to 0} \frac{1}{h^{-\alpha}} \sum_{m=0}^{\infty} (-1)^{m}
\frac{\Gamma(-\alpha + 1)}{\Gamma(m + 1)\Gamma(-\alpha - m + 1)} f(x - mh)
\end{aligned}
$$

[負の整数について一般化された二項係数](https://en.wikipedia.org/wiki/Binomial_coefficient#Generalization_to_negative_integers)を使って式を変形します。

 負の整数について一般化された二項係数 $\binom{-n}{m}$ を展開します。

$$
\begin{aligned}
\binom{-n}{m}
&= \frac{-n(-n - 1)(-n - 2)(-n - 3) \ldots (-n - m + 1)}{m!}\\
&= (-1)^{m} \frac{n(n + 1)(n + 2)(n + 3) \ldots (n + m - 1)}{m!}\\
&= (-1)^{m} \frac{(n + m - 1)!}{m!(n - 1)!}\\
&= (-1)^{m} \frac{\Gamma(n + m)}{\Gamma(m + 1)\Gamma(n)}
\end{aligned}
$$

これで $\alpha$ 階積分が定義できます。

$$
d^{-\alpha} f(x)
= \lim_{h \to 0} \frac{1}{h^{-\alpha}} \sum_{m=0}^{\infty}
\frac{\Gamma(\alpha + m)}{\Gamma(m + 1)\Gamma(\alpha)} f(x - mh)
\,, \qquad \alpha \in \Bbb{R}
$$

## 分数階微積分の実装
Grünwald–Letnikovの分数階微分と分数階積分を再掲します。

$$
\begin{aligned}
d^{\alpha} f(x)
&= \lim_{h \to 0} \frac{1}{h^{\alpha}} \sum_{m=0}^{\infty} (-1)^{m}
\frac{\Gamma(\alpha + 1)}{\Gamma(m + 1)\Gamma(\alpha - m + 1)} f(x - mh)\\

d^{-\alpha} f(x)
&= \lim_{h \to 0} \frac{1}{h^{-\alpha}} \sum_{m=0}^{\infty}
\frac{\Gamma(\alpha + m)}{\Gamma(m + 1)\Gamma(\alpha)} f(x - mh)
\end{aligned}
$$

実装します。

```javascript
// ガンマ関数の計算で定義した関数 gamma を使用。

function fractionalDerivative(func, x, alpha, h) {
  var sum = 0
  var sign = 1
  for (var m = 0; m < MAX_SUM; ++m) {
    sum += sign / gamma(m + 1) / gamma(alpha - m + 1) * func(x - m * h)
    sign *= -1
  }
  return sum * gamma(alpha + 1) / h ** alpha
}

function fractionalIntegral(func, x, alpha, h) {
  alpha = -alpha
  var sum = 0
  for (var m = 0; m < MAX_SUM; ++m) {
    sum += gamma(alpha + m) / gamma(m + 1) * func(x - m * h)
  }
  return sum / gamma(alpha) * h ** alpha
}

// 計算する関数とパラメータ。
// alpha が正の値なら微分、負の値なら積分。
var MAX_SUM = 128
var exampleFunc = (x) => { return x * Math.sin(x) }
var alpha = 0.5
var h = 0.01

var maxI = 1024 - 1
var maxX = 10
var result = []
for (var i = 0; i <= maxI; ++i) {
  var x = maxX * i / maxI
  var deriv = fractionalDerivative(exampleFunc, x, alpha, h)
  var integ = fractionalIntegral(exampleFunc, x, alpha, h)
  result.push([deriv, integ])
}
console.log(result)
```

## デモ
<script src="demo/grunwald_letnikov/fft.js"></script>
<script src="demo/grunwald_letnikov/vec2.js"></script>
<script src="demo/grunwald_letnikov/canvas.js"></script>
<script>
function gamma(z) {
  var p = [
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  ]

  if (z < 0.5) {
    return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z))
  }

  z -= 1
  var a = 0.99999999999980993
  for (var i = 0; i < p.length; ++i) {
    a += p[i] / (z + i + 1)
  }
  var t = z + p.length - 0.5
  return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * a
}

class GrunwaldLetnikov {
  constructor(func, alpha, delta, maxX = 1) {
    this.func = func
    this.alpha = alpha
    this.delta = delta
    this.maxIter = 128

    this.maxX = maxX
  }

  fractionalDerivative(x) {
    var sum = 0
    var sign = 1
    for (var m = 0; m < this.maxIter; ++m) {
      sum += sign / gamma(m + 1) / gamma(this.alpha - m + 1)
        * this.func(x - m * this.delta)
      sign *= -1
    }
    return sum * gamma(this.alpha + 1) / this.delta ** this.alpha
  }

  fractionalIntegral(x) {
    var alpha = -this.alpha
    var sum = 0
    for (var m = 0; m < this.maxIter; ++m) {
      sum += gamma(alpha + m) / gamma(m + 1) * this.func(x - m * this.delta)
    }
    return sum / gamma(alpha) * this.delta ** alpha
  }

  fractionalDifferintegral(x) {
    if (this.alpha > 0) {
      return this.fractionalDerivative(x)
    }
    else if (this.alpha < 0) {
      return this.fractionalIntegral(x)
    }
    return func(x)
  }

  getDerivative(canvas, func) {
    var line = new Array(canvas.width).fill(0)
    var min = Number.MAX_VALUE
    var max = -Number.MAX_VALUE
    var denom = canvas.width - 1
    for (var x = 0; x < canvas.width; ++x) {
      line[x] = new Vec2(x, func(2 * this.maxX * x / denom - this.maxX))
      if (max < line[x].y) {
        max = line[x].y
      }
      if (min > line[x].y) {
        min = line[x].y
      }
    }

    min = Number.isFinite(min) ? min : Number.MAX_VALUE
    max = Number.isFinite(max) ? max : -Number.MAX_VALUE
    return [line, min, max]
  }

  drawText(canvas, min, max) {
    var fontSize = 12
    var fontSizeHalf = fontSize / 2
    var lineSpace = fontSize + 2
    canvas.context.font = `${fontSize}px monospace`

    var texts = [
      `x: [0, ${this.maxX}]`,
      `y: [${min}, ${max}]`
    ]

    var maxTextWidth = 0
    for (var i = 0; i < texts.length; ++i) {
      var textWidth = canvas.context.measureText(texts[i]).width
      if (maxTextWidth < textWidth) {
        maxTextWidth = textWidth
      }
    }

    canvas.context.fillStyle = "#ffffffb0"
    canvas.context.fillRect(0, 0, maxTextWidth + fontSizeHalf, 2.2 * lineSpace)

    canvas.context.fillStyle = "#303030"
    for (var i = 0; i < texts.length; ++i) {
      canvas.context.fillText(texts[i], fontSizeHalf, (i + 1) * lineSpace)
    }
  }

  drawPlot(canvas) {
    var [lineD, minD, maxD]
      = this.getDerivative(canvas, (x) => this.fractionalDerivative(x))
    var [lineI, minI, maxI]
      = this.getDerivative(canvas, (x) => this.fractionalIntegral(x))
    var min = minD < minI ? minD : minI
    var max = maxD > maxI ? maxD : maxI
    var diff = max - min
    for (var i = 0; i < canvas.width; ++i) {
      lineD[i].y = canvas.height - canvas.height * (lineD[i].y - min) / diff
      lineI[i].y = canvas.height - canvas.height * (lineI[i].y - min) / diff
    }

    var yZero = canvas.height + min / diff * canvas.height
    canvas.context.lineWidth = 1
    canvas.context.setLineDash([])
    canvas.context.strokeStyle = "#808080"
    canvas.drawLine(new Vec2(0, yZero), new Vec2(canvas.width, yZero))
    var halfWidth = canvas.width / 2
    canvas.drawLine(new Vec2(halfWidth, 0), new Vec2(halfWidth, canvas.height))

    canvas.context.font = `18px monospace`
    canvas.context.fillStyle = "#303030"
    canvas.context.fillText(`0`, halfWidth, yZero - 2)

    canvas.context.lineJoin = "round"
    canvas.context.lineWidth = 3

    canvas.context.setLineDash([])
    canvas.context.strokeStyle = "#303030"
    canvas.drawPath(lineD)

    canvas.context.setLineDash([1])
    canvas.context.strokeStyle = "#ff0000"
    canvas.drawPath(lineI)

    return [min, max]
  }

  draw(canvas) {
    canvas.clearWhite()

    // Plotting.
    var [min, max] = this.drawPlot(canvas)
    this.drawText(canvas, min, max)
  }
}

// Fourier transform based derivative computation.
// This class uses FFT.js by indutny.
// https://github.com/indutny/fft.js
class FtDerivative {
  constructor(func = (x) => Math.sin(2 * Math.PI * x), alpha = 1) {
    this.alpha = alpha

    this.data = new Array(ftDerivativeCanvas.width)
    var dataLast = this.data.length - 1
    for (var i = 0; i < this.data.length; ++i) {
      this.data[i] = func(i / dataLast)
    }
    this.fft = new FFT(this.data.length)
    this.inputSpectrum = this.fft.createComplexArray()
    this.fft.realTransform(this.inputSpectrum, this.data)
    this.specLast = this.inputSpectrum.length - 1

    this.outputSpectrum = new Array(this.inputSpectrum.length).fill(0)
    this.outputComplex = new Array(this.inputSpectrum.length).fill(0)
    this.output = new Array(this.data.length)
    this.setMinMax(this.output)
  }

  compute() {
    this.outputSpectrum[0] = this.inputSpectrum[0]
    this.outputSpectrum[1] = this.inputSpectrum[1]
    for (var curr = 2; curr < this.inputSpectrum.length; curr += 2) {
      var next = curr + 1

      // (j * omega ** alpha) * F(omega)
      var omega = Math.PI * curr / this.specLast
      var real = this.alpha * Math.log(Math.abs(omega))
      var imag = this.alpha * Math.PI / 2

      var expReal = Math.exp(real)
      real = expReal * Math.cos(imag)
      imag = expReal * Math.sin(imag)

      var specReal = this.inputSpectrum[curr]
      var specImag = this.inputSpectrum[next]
      var realTemp = real
      this.outputSpectrum[curr] = specReal * real - specImag * imag
      this.outputSpectrum[next] = specReal * imag + specImag * realTemp
    }
    this.fft.inverseTransform(this.outputComplex, this.outputSpectrum)
    this.output = this.fft.fromComplexArray(this.outputComplex)
    this.setMinMax(this.output)
  }

  setMinMax(data) {
    this.min = Number.MAX_VALUE
    this.max = -Number.MAX_VALUE
    for (var i = 0; i < data.length; ++i) {
      if (this.max < data[i]) {
        this.max = data[i]
      }
      if (this.min > data[i]) {
        this.min = data[i]
      }
    }
    this.diff = this.max - this.min
  }

  draw(canvas) {
    this.compute()

    var line = new Array(canvas.width)
    for (var x = 0; x < canvas.width; ++x) {
      line[x] = new Vec2(
        x,
        (this.output[x] - this.min) / this.diff * canvas.height
      )
    }
    canvas.clearWhite()
    canvas.drawPath(line)
  }
}

function refresh() {
  glDerivative.alpha = alpha
  glDerivative.func = func
  glDerivative.draw(glDerivativeCanvas)

  // ftDerivative.alpha = alpha
  // ftDerivative.func = func
  // ftDerivative.draw(ftDerivativeCanvas)
}

var func = (x) => Math.sin(2 * Math.PI * x)
var alpha = 1

var glDerivativeCanvas = new Canvas(document.body, 640, 320)
var glDerivative = new GrunwaldLetnikov(
  func, alpha, 1 / glDerivativeCanvas.width, 1)

// var ftDerivativeCanvas = new Canvas(document.body, 512, 256)
// var ftDerivative = new FtDerivative(func, alpha)

var numberAlpha = new NumberInput(
  document.body, "alpha", alpha, -10, 10, 0.01,
  {
    event: "input", func: (value) => { alpha = value; refresh() },
  }
)
var numberH = new NumberInput(
  document.body, "h", glDerivative.delta, 0.001, 1, 0.001,
  {
    event: "input",
    func: (value) => {
      glDerivative.delta = value
      glDerivative.draw(glDerivativeCanvas)
    },
  }
)

refresh()
</script>

黒のプロットが `fractionalDerivative` 、赤が `fractionalIntegral` です。

ブラウザの開発者コンソールから関数 `func` を上書きすることで描画する関数を変更できます。 `func` の第一引数には横軸の座標が入力されます。コンソールから再描画するには関数 `refresh` を呼び出してください。

```javascript
func = (x) => x**5
refresh()
```

## 問題点
$0 < \alpha < 1$ の範囲を超えると正しく計算できていないように見えます。この問題は $\alpha$ を整数部と小数部に分けて計算すれば解決できるようです。例えば2.3階微分を計算するときは、何か別の方法で2階微分してから、Grünwald–Letnikovで0.3階微分します。

## 参考文献
- [Fractional calculus - Wikipedia](https://en.wikipedia.org/wiki/Fractional_calculus)
- [Fractional Calculus: History, Definitions and Applications for the Engineer](https://web.archive.org/web/20051029113800/http://www.nd.edu/~msen/Teaching/UnderRes/FracCalc.pdf)
- [Fractional Calculus - Xuru's Website](http://www.xuru.org/fc/TOC.asp)
- [Gamma function - Rosetta Code](https://rosettacode.org/wiki/Gamma_function)
- [Fractional order derivative of a function & fractional numbers' factorial. - YouTube](https://www.youtube.com/watch?v=7PCQmxlX9mU)
- [The Grünwald–Letnikov method for fractional differential equations - ScienceDirect](https://www.sciencedirect.com/science/article/pii/S0898122111002173)
