class Burgers1DInviscid {
  // 波はインデックス x が大きくなる方向に向かって進む。
  constructor(length, dx, dt) {
    this.length = length
    this.reset()

    this.dx = dx
    this.dt = dt

    this.refreshConstants()
    this.pick(0, 0)
  }

  refreshConstants() {
    this.C0 = this.dt / this.dx / 2
  }

  step() {
    this.wave.unshift(this.wave.pop())

    // Left boundary.
    this.wave[0][0] = 0

    for (var x = 1; x < this.length; ++x) {
      var left = x - 1
      this.wave[0][x] = this.wave[1][x] - this.C0 * (
        this.wave[1][x] * this.wave[1][x]
        - this.wave[1][left] * this.wave[1][left]
      )
    }

    if (this.pickY != 0) {
      this.wave[0][this.pickX] = this.pickY
    }
  }

  reset() {
    this.wave = new Array(2)
    this.wave[0] = new Array(this.length).fill(0)
    this.wave[1] = new Array(this.length).fill(0)
  }

  pick(x, y) {
    x = Math.min(Math.max(x, 0.0), 1.0)
    this.pickX = Math.floor((this.length - 1) * x)
    this.pickY = y
  }
}

class Burgers1DInviscidGodunov {
  // pick に入力できる y の最大値は dx / dt 。
  constructor(length, dx, dt) {
    this.length = length
    this.reset()

    this.dx = dx
    this.dt = dt

    this.refreshConstants()
    this.pick(0, 0)
  }

  refreshConstants() {
    this.C0 = this.dt / this.dx / 2
  }

  uStar(u, v) {
    if (u >= v) return (u + v) / 2 > 0 ? u : v
    else if (u > 0) return u
    else if (v < 0) return v
    return 0
  }

  step() {
    this.wave.unshift(this.wave.pop())

    var last = this.length - 1

    // Left boundary.
    this.wave[0][0] = 0
    this.wave[0][last] = 0

    for (var x = 1; x < last; ++x) {
      var uStarL = this.uStar(this.wave[1][x], this.wave[1][x + 1])
      var uStarR = this.uStar(this.wave[1][x - 1], this.wave[1][x])
      this.wave[0][x] = this.wave[1][x] - this.C0 * (uStarL * uStarL - uStarR * uStarR)
    }

    if (this.pickY != 0) {
      this.wave[0][this.pickX] = this.pickY
    }
  }

  reset() {
    this.wave = new Array(2)
    this.wave[0] = new Array(this.length).fill(0)
    this.wave[1] = new Array(this.length).fill(0)
  }

  pick(x, y) {
    x = Math.min(Math.max(x, 0.0), 1.0)
    this.pickX = Math.floor((this.length - 1) * x)
    this.pickY = y
  }
}

class Burgers1DViscid {
  // 波はインデックス x が大きくなる方向に向かって進む。
  constructor(length, dx, dt, epsilon = 0.01) {
    this.length = length
    this.reset()

    this.dx = dx
    this.dt = dt
    this.epsilon = epsilon

    this.refreshConstants()
    this.pick(0, 0)
  }

  refreshConstants() {
    this.C1 = this.dt / this.dx
    this.C0 = this.C1 * this.epsilon / this.dx / 2
  }

  step() {
    this.wave.unshift(this.wave.pop())

    var last = this.length - 1

    // Boundary.
    this.wave[0][0] = 0
    this.wave[0][last] = 0

    for (var x = 1; x < last; ++x) {
      var left = x - 1
      var right = x + 1
      this.wave[0][x] = this.wave[1][x]
        + this.C0 * (
          this.wave[1][right] - 2 * this.wave[1][x] + this.wave[1][left])
        + this.C1 * (
          this.wave[1][left] - this.wave[1][right])
    }

    if (this.pickY != 0) {
      this.wave[0][this.pickX] = this.pickY
    }
  }

  reset() {
    this.wave = new Array(2)
    this.wave[0] = new Array(this.length).fill(0)
    this.wave[1] = new Array(this.length).fill(0)
  }

  pick(x, y) {
    x = Math.min(Math.max(x, 0.0), 1.0)
    this.pickX = Math.floor((this.length - 1) * x)
    this.pickY = y
  }
}

/* UI */
function onMouseDownCanvas(event) {
  mousedown = true
  var rect = event.target.getBoundingClientRect()
  var x = event.clientX - rect.left
  var y = event.clientY - rect.top

  wave.pick(x / cv.width, cv.height / 2 - y)
}

function onMouseMoveCanvas(event) {
  if (!mousedown) return
  var rect = event.target.getBoundingClientRect()
  var x = event.clientX - rect.left
  var y = event.clientY - rect.top

  wave.pick(x / cv.width, cv.height / 2 - y)
}

function onMouseUpCanvas(event) {
  mousedown = false
  wave.pick(0.5, 0)
}

function onMouseOutCanvas(event) {
  onMouseUpCanvas(event)
}

function updateCanvas() {
  cv.clearWhite()

  xCenter = cv.width / 2
  yCenter = cv.height / 2

  cv.context.save()
  cv.context.translate(0, yCenter)
  cv.context.scale(1, -1)

  position = Array.isArray(wave.wave[0]) ? wave.wave[0] : wave.wave

  cv.context.lineWidth = 10
  cv.context.strokeStyle = "#4488ff"
  cv.context.fillStyle = "#4488ff"
  cv.context.lineJoin = "round"
  cv.context.beginPath()
  cv.context.moveTo(0, position[0])
  for (let x = 1; x < width; ++x) {
    cv.context.lineTo(x, position[x])
  }
  cv.context.lineTo(width, position[width - 1])
  cv.context.lineTo(width, -cv.height)
  cv.context.lineTo(0, -cv.height)
  cv.context.closePath()
  cv.context.fill()

  cv.context.restore()

  cv.context.lineWidth = 1
  cv.context.strokeStyle = "#000000"
  cv.context.beginPath()
  cv.context.moveTo(0, yCenter)
  cv.context.lineTo(cv.width, yCenter)
  cv.context.stroke()
}

class Lowpass {
  constructor(samplerate, cutoff) {
    this.samplerate = samplerate
    this.cutoff = cutoff
    this.s = 0
  }

  set cutoff(freq) {
    var omega_c = Math.tan(Math.PI * freq / this.samplerate)
    this.g = omega_c / (1 + omega_c)
  }

  process(x) {
    var v = (x - this.s) * this.g
    var y = v + this.s
    this.s = y + v
    return y
  }
}

function stepWave() {
  for (var i = 0; i < 16; ++i) {
    wave.step()
  }
}

function animate() {
  updateCanvas()
  stepWave()
  requestAnimationFrame(animate)
}

/* main */
var stepPerFrame = 4

var width = 512
var height = 512

var dx = 1
var dt = 1 / 512

// var wave = new Burgers1DInviscid(width, dx, dt)
var wave = new Burgers1DInviscidGodunov(width, dx, dt)
// var wave = new Burgers1DViscid(width, dx, dt, 1e-5)

var divWave = new Div(document.body, "divWave")
var cv = new Canvas(divWave.element, width, height)
var buttonReset = new Button(divWave.element,
  "Reset",
  () => { wave.reset() })

var mousedown = false
cv.element.addEventListener("mousedown", onMouseDownCanvas, false)
cv.element.addEventListener("mousemove", onMouseMoveCanvas, false)
cv.element.addEventListener("mouseup", onMouseUpCanvas, false)
cv.element.addEventListener("mouseout", onMouseOutCanvas, false)

animate()
