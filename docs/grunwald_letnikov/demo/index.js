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

var func = (x) => Math.sin(2 * Math.PI * x ** 2)
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
