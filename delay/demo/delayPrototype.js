class LagrangeNInterp {
  constructor(order) {
    this.N = order
    this.fix = (this.N - 1) / 2
    this.xd = new Array(this.N).fill(0)
    this.diff = new Array(this.N).fill(0)
  }

  reset() {
    this.xd.fill(0)
    this.diff.fill(0)
  }

  push(input) {
    this.diff[0] = input - this.xd[0]
    this.xd[0] = input
    for (var i = 1; i < this.N; ++i) {
      this.diff[i] = this.diff[i - 1] - this.xd[i]
      this.xd[i] = this.diff[i - 1]
    }
  }

  at(fraction) {
    var delta = fraction + this.fix
    var sig = 0

    var i = this.N
    while (i > 0) {
      var next = i - 1
      sig = (next - delta) / i * (this.diff[next] + sig)
      i = next
    }

    return sig + this.xd[0]
  }

  process(input, fraction) {
    this.push(input)
    return this.at(fraction)
  }
}

class Delay {
  // readOrder と writeOrder は奇数。
  constructor(sampleRate, time, overSample = 2, readOrder = 1, writeOrder = 1) {
    this.overSample = overSample
    this.sampleRate = this.overSample * sampleRate

    // 補間で追加されるディレイの補正値。
    this.fix = (readOrder - 1 + this.overSample * (writeOrder - 1)) / 2

    this.buf = new Array(Math.ceil(time * this.sampleRate) + 1).fill(0)

    this.rInterp = new LagrangeNInterp(readOrder)
    this.wInterp = new LagrangeNInterp(writeOrder)

    this.wptr = 0
    this.time = time
  }

  set time(value) {
    var delayTime = Math.max(
      this.fix, Math.min(this.sampleRate * value, this.buf.length))
    this.rFraction = delayTime % 1
    this.rptr = Math.ceil(this.wptr - delayTime) + this.fix

    // mod(rptr, buf.length). 正の余りを計算。
    while (this.rptr < 0) this.rptr += this.buf.length
    this.rptr %= this.buf.length
  }

  reset() {
    this.buf.fill(0)
    this.rInterp.reset()
    this.wInterp.reset()
  }

  process(input) {
    this.wInterp.push(input)

    for (var i = 1; i <= this.overSample; ++i) {
      var wFraction = (this.overSample - i) / this.overSample

      this.buf[this.wptr] = this.wInterp.at(wFraction)
      this.wptr = (this.wptr + 1) % this.buf.length

      this.rInterp.push(this.buf[this.rptr])
      this.rptr = (this.rptr + 1) % this.buf.length
    }

    return this.rInterp.at(this.rFraction)
  }
}

/// テスト。
function generateSin(sampleRate, length, frequency) {
  var wave = new Array(sampleRate * length).fill(0)
  var sin_omega_fs = 2 * Math.PI * frequency / sampleRate
  for (var i = 0; i < wave.length; ++i) {
    wave[i] = Math.sin(i * sin_omega_fs)
  }
  return wave
}

function generateNoise(sampleRate, length) {
  var wave = new Array(sampleRate * length).fill(0)
  for (var i = 0; i < wave.length; ++i) {
    wave[i] = Math.random()
  }
  return wave
}

function renderDelay(delay, source, lfo, modBase, modAmount) {
  var dest = new Array(source.length).fill(0)
  var fraction = new Array(source.length).fill(0)
  for (var i = 0; i < dest.length; ++i) {
    dest[i] = delay.process(source[i])
    fraction[i] = delay.fraction
    // delay.time = modBase + modAmount * Math.exp(lfo[i] - 1)
    delay.time = modBase + modAmount * lfo[i]
  }
  return [dest, fraction]
}

function diffInt(sampleRate, sig, modBase, modAmount) {
  var diff = new Array(sig.length).fill(0)
  var t1 = modBase
  for (var i = 0; i < sig.length; ++i) {
    var t0 = modBase + modAmount * (1 - Math.exp(sig[i] - 1))
    diff[i]
      = Math.floor(sampleRate * t0)
      - Math.floor(sampleRate * t1)
    t1 = t0
  }
  return diff
}

function testModulation() {
  var sampleRate = 44100
  var dur = 2 // 秒。
  var modBase = 0.001
  var modAmount = 1.9

  var result = {}

  // var source = generateNoise(sampleRate, dur)
  var source = generateSin(sampleRate, dur, 220)
  var lfo = generateSin(sampleRate, dur, 2)

  var denom = lfo.length - 1
  for (var i = 0; i < lfo.length; ++i) {
    lfo[i] = 0 * lfo[i] + 1.00 * i / denom
  }

  result["source"] = { sampleRate, wave: source }
  // result["lfo"] = { sampleRate, wave: lfo }
  // result["lfoDiff"] = { sampleRate, wave: diffInt(sampleRate, lfo, modBase, modAmount) }

  for (var order = 1; order <= 9; order += 2) {
    var overSample = 2//order == 1 ? 2 : order
    var rOrder = order
    var wOrder = order
    var delay = new Delay(sampleRate, dur, overSample, order, wOrder)

    var [dest, fraction] = renderDelay(delay, source, lfo, modBase, modAmount)

    var overSampleStr = String(overSample).padStart(2, "0")
    result[`${overSampleStr}x_r${rOrder}_w${wOrder}`]
      = { sampleRate, wave: dest }
    // result[`${overSampleStr}x_r${rOrder}_w${wOrder}_fraction`]
    //   = { sampleRate, wave: fraction }
  }

  console.log(result)
}

function renderDelayIR(sampleRate, delay, nTest = 11) {
  // 16 サンプルのディレイ。
  var impulse = new Array(32).fill(0)
  impulse[0] = 1

  var dest = new Array(nTest)
  for (var n = 0; n < dest.length; ++n) {
    dest[n] = new Array(impulse.length)
    delay.time = (16 + n / (nTest - 1)) / sampleRate
    delay.reset()
    for (var i = 0; i < impulse.length; ++i) {
      dest[n][i] = delay.process(impulse[i])
    }
  }
  return dest
}

function testDelayIR() {
  var sampleRate = 44100
  var maxTime = 0.1

  var result = {}

  var overSample = 9
  for (var order = 1; order <= 9; order += 2) {
    var rOrder = order
    var wOrder = order

    var delay = new Delay(sampleRate, maxTime, overSample, rOrder, wOrder)
    var overSampleStr = String(overSample).padStart(2, "0")
    result[`${overSampleStr}x_r${rOrder}_w${wOrder}`] = renderDelayIR(sampleRate, delay)
  }

  console.log(result)
}

// testModulation()
// testDelayIR()
