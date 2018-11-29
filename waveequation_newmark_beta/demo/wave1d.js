/* math functions */
class Solver {
  // Ax = b の A が固定で 0 が多く含まれているなら速い。
  constructor(a = []) {
    this.prepareA(a)
    this.tolerance = 1e-9
    this.maxIteration = 1024
  }

  prepareA(a) {
    this.aReduced = new Array(a.length)
    for (var i = 0; i < a.length; ++i) {
      this.aReduced[i] = []
      for (var j = 0; j < a[i].length; ++j) {
        if (a[i][j] === 0 || i === j) {
          continue
        }
        this.aReduced[i].push([a[i][j], j])
      }
    }

    this.aDiag = new Array(a.length)
    for (var i = 0; i < a.length; ++i) {
      this.aDiag[i] = a[i][i]
    }

    this.x = new Array(a.length).fill(0)
    this.x_prev = new Array(a.length).fill(0)
  }

  solve(b) {
    this.x.fill(0)
    for (var iter = 0; iter < this.maxIteration; ++iter) {
      this.x_prev = [...this.x]
      for (var i = 0; i < this.x.length; ++i) {
        var sum = b[i]
        for (var j = 0; j < this.aReduced[i].length; ++j) {
          sum -= this.aReduced[i][j][0] * this.x[this.aReduced[i][j][1]]
        }
        this.x[i] = sum / this.aDiag[i]
      }
      if (this.isInTolerance(this.x)) {
        break
      }
    }
    return this.x
  }

  solveBX(b, x) {
    x.fill(0)
    for (var iter = 0; iter < this.maxIteration; ++iter) {
      this.x_prev = [...x]
      for (var i = 0; i < x.length; ++i) {
        var sum = b[i]
        for (var j = 0; j < this.aReduced[i].length; ++j) {
          sum -= this.aReduced[i][j][0] * x[this.aReduced[i][j][1]]
        }
        x[i] = sum / this.aDiag[i]
      }
      if (this.isInTolerance(x)) {
        break
      }
    }
    return x
  }

  isInTolerance(vecX) {
    for (var i = 0; i < vecX.length; ++i) {
      if (Math.abs(vecX[i] - this.x_prev[i]) > this.tolerance) return false
    }
    return true
  }

  solveAB(a, b) {
    this.prepareA(a)
    return this.solve(b)
  }
}

/* wave */
class Wave1DNewmarkBeta {
  constructor(length, c, dx, dt, attenuation) {
    this.length = length
    this.reset()

    this.c = c
    this.dx = dx
    this.dt = dt
    this.attenuation = attenuation

    this.beta = 1 / 4
    this.vectorB = new Array(this.length).fill(0)

    this.boundary = 2 // 1: 固定端, 2: 自由端。

    this.refreshConstants()
    this.pick(0, 0)
  }

  refreshConstants() {
    var dt2 = this.dt * this.dt

    this.C0 = (this.dx / this.c) ** 2 + 2 * dt2 * this.beta
    this.C1 = - dt2 * this.beta
    this.C2 = dt2 / 2 + this.C1
    this.C3 = this.dt / 2

    this.L = this.boundary * this.C1
    this.R = this.boundary * this.C1

    this.initMatrix()
  }

  initMatrix() {
    var mat = []
    for (var i = 0; i < this.length; ++i) {
      mat.push(new Array(this.length).fill(0))
    }

    mat[0][0] = this.C0
    mat[0][1] = this.L

    var last = this.length - 1
    for (var i = 1; i < last; ++i) {
      mat[i][i - 1] = this.C1
      mat[i][i] = this.C0
      mat[i][i + 1] = this.C1
    }

    mat[last][last - 1] = this.R
    mat[last][last] = this.C0

    if (this.solver === undefined) {
      this.solver = new Solver(mat)
      this.solver.tolerance = 1e-9
      this.solver.maxIteration = 1024
    }
    else {
      this.solver.prepareA(mat)
    }
  }

  step() {
    this.acceleration.unshift(this.acceleration.pop())

    var last = this.length - 1

    // 左端。
    var right = 1
    this.vectorB[0] = (this.boundary * this.wave[right] - 2 * this.wave[0])
      + this.dt * (this.boundary * this.velocity[right] - 2 * this.velocity[0])
      + this.C2 * (this.boundary * this.acceleration[1][right]
        - 2 * this.acceleration[1][0])

    // 右端。
    var left = last - 1
    this.vectorB[last] = (this.boundary * this.wave[left] - 2 * this.wave[last])
      + this.dt * (this.boundary * this.velocity[left] - 2 * this.velocity[last])
      + this.C2 * (this.boundary * this.acceleration[1][left]
        - 2 * this.acceleration[1][last])

    for (var x = 1; x < last; ++x) {
      left = x - 1
      right = x + 1
      this.vectorB[x] = (this.wave[left] - 2 * this.wave[x] + this.wave[right])
        + this.dt * (this.velocity[left]
          - 2 * this.velocity[x] + this.velocity[right])
        + this.C2 * (this.acceleration[1][left]
          - 2 * this.acceleration[1][x] + this.acceleration[1][right])
    }

    // c * dt / dx が非常に大きい場合にソルバが収束しない。
    // solver.maxIteration の値を上げれば計算コストと引きかえに発散を抑えることができる。
    this.solver.solveBX(this.vectorB, this.acceleration[0])

    for (var x = 0; x < this.length; ++x) {
      this.wave[x] += this.dt * this.velocity[x]
        + this.C2 * this.acceleration[1][x]
        - this.C1 * this.acceleration[0][x]
      this.velocity[x] += this.C3
        * (this.acceleration[1][x] + this.acceleration[0][x])
      this.wave[x] *= this.attenuation
    }

    if (this.pickY != 0) {
      this.wave[this.pickX] = this.pickY
    }
  }

  reset() {
    this.wave = new Array(this.length).fill(0)
    this.velocity = new Array(this.length).fill(0)

    this.acceleration = []
    for (var i = 0; i < 2; ++i) {
      this.acceleration.push(new Array(this.length).fill(0))
    }
  }

  pick(x, y) {
    x = Math.min(Math.max(x, 0.0), 1.0)
    this.pickX = Math.floor((this.length - 1) * x)
    this.pickY = y
  }
}

class DampedWave1DNewmarkBeta {
  constructor(length, c, dx, dt, a, k) {
    this.length = length
    this.reset()

    this.c = c
    this.dt = dt
    this.dx = dx
    this.a = a
    this.k = k

    this.beta = 1 / 4
    this.vectorB = new Array(this.length).fill(0)

    this.boundary = 1 // 1: 固定端, 2: 自由端。

    this.refreshConstants()
    this.pick(0, 0)
  }

  refreshConstants() {
    var dt2 = this.dt * this.dt

    this.C2 = this.c * this.c / dt2
    this.C3 = dt2 * (0.5 - this.beta)
    this.C7 = this.dt / 2
    this.C8 = dt2 * this.beta
    this.C1 = - this.C2 * this.C8
    this.C6 = this.k + 2 * this.C2
    this.C0 = 1 + this.a * this.C7 + this.C6 * this.C8
    this.C4 = this.C3 * this.C6 + this.a * this.C7
    this.C5 = this.a + this.dt * this.C6

    this.L = this.boundary * this.C1
    this.R = this.boundary * this.C1

    this.initMatrix()
  }

  initMatrix() {
    var mat = []
    for (var i = 0; i < this.length; ++i) {
      mat.push(new Array(this.length).fill(0))
    }

    mat[0][0] = this.C0
    mat[0][1] = this.L

    var last = this.length - 1
    for (var i = 1; i < last; ++i) {
      mat[i][i - 1] = this.C1
      mat[i][i] = this.C0
      mat[i][i + 1] = this.C1
    }

    mat[last][last - 1] = this.R
    mat[last][last] = this.C0

    if (this.solver === undefined) {
      this.solver = new Solver(mat)
      this.solver.tolerance = 1e-9
      this.solver.maxIteration = 1024
    }
    else {
      this.solver.prepareA(mat)
    }
  }

  step() {
    this.acceleration.unshift(this.acceleration.pop())

    var last = this.length - 1

    // 左端。
    var right = 1
    this.vectorB[0]
      = this.boundary * this.C2 * (
        this.wave[right]
        + this.dt * this.velocity[right]
        + this.C3 * this.acceleration[1][right]
      )
      - this.C4 * this.acceleration[1][0]
      - this.C5 * this.velocity[0]
      - this.C6 * this.wave[0]

    // 右端。
    var left = last - 1
    this.vectorB[last]
      = this.boundary * this.C2 * (
        this.wave[left]
        + this.dt * this.velocity[left]
        + this.C3 * this.acceleration[1][left]
      )
      - this.C4 * this.acceleration[1][last]
      - this.C5 * this.velocity[last]
      - this.C6 * this.wave[last]

    for (var x = 1; x < last; ++x) {
      left = x - 1
      right = x + 1
      this.vectorB[x]
        = this.C2 * (
          (this.wave[left] + this.wave[right])
          + this.dt * (this.velocity[left] + this.velocity[right])
          + this.C3 * (this.acceleration[1][left] + this.acceleration[1][right])
        )
        - this.C4 * this.acceleration[1][x]
        - this.C5 * this.velocity[x]
        - this.C6 * this.wave[x]
    }

    // c * dt / dx が非常に大きい場合にソルバが収束しない。
    // solver.maxIteration の値を上げれば計算コストと引きかえに発散を抑えることができる。
    this.solver.solveBX(this.vectorB, this.acceleration[0])

    for (var x = 0; x < this.length; ++x) {
      this.wave[x] += this.dt * this.velocity[x]
        + this.C3 * this.acceleration[1][x]
        + this.C8 * this.acceleration[0][x]
      this.velocity[x] += this.C7
        * (this.acceleration[1][x] + this.acceleration[0][x])
    }

    if (this.pickY != 0) {
      this.wave[this.pickX] = this.pickY
    }
  }

  reset() {
    this.wave = new Array(this.length).fill(0)
    this.velocity = new Array(this.length).fill(0)

    this.acceleration = []
    for (var i = 0; i < 2; ++i) {
      this.acceleration.push(new Array(this.length).fill(0))
    }
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

  wave.pick(x / cv.width, y - cv.height / 2)
}

function onMouseMoveCanvas(event) {
  if (!mousedown)
    return
  var rect = event.target.getBoundingClientRect()
  var x = event.clientX - rect.left
  var y = event.clientY - rect.top

  wave.pick(x / cv.width, y - cv.height / 2)
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
  cv.context.lineTo(cv.width, cv.height)
  cv.context.lineTo(0, cv.height)
  cv.context.closePath()
  cv.context.fill()

  cv.context.lineWidth = 0.1
  cv.context.strokeStyle = "#404040"
  cv.context.beginPath()
  cv.context.moveTo(0, yCenter)
  cv.context.lineTo(cv.width, yCenter)
  cv.context.stroke()

  cv.context.restore()
}

function animate() {
  updateCanvas()
  wave.step()
  requestAnimationFrame(animate)
}

/* main */
var stepPerFrame = 4

var width = 512
var height = 512

// var wave = new Wave1DNewmarkBeta(width, 4, 0.1, 1 / 60, 1)
var wave = new DampedWave1DNewmarkBeta(width, 4, 0.1, 1 / 60, 1, 1)

var divWave = new Div(document.body, "divWave")
var cv = new Canvas(divWave.element, width, height)
var buttonReset = new Button(divWave.element,
  "Reset",
  () => { wave.reset() })
var numberC = new NumberInput(divWave.element,
  "c", wave.c, 0.001, 128, 0.001,
  (value) => { wave.c = value; wave.refreshConstants() })
var numberDx = new NumberInput(divWave.element,
  "Δ_x", wave.dx, 0.01, 1, 0.01,
  (value) => { wave.dx = value; wave.refreshConstants() })
var numberA = new NumberInput(divWave.element,
  "a", wave.a, 0, 128, 0.001,
  (value) => { wave.a = value; wave.refreshConstants() })
var numberK = new NumberInput(divWave.element,
  "k", wave.k, 0, 128, 0.001,
  (value) => { wave.k = value; wave.refreshConstants() })
var numberBeta = new NumberInput(divWave.element,
  "β", wave.beta, 0, 2, 0.001,
  (value) => { wave.beta = value; wave.refreshConstants() })
var pullDownMenuBoundaryCondition = new PullDownMenu(divWave.element,
  "Boundary Condition",
  (value) => {
    if (value === "Free") {
      wave.boundary = 2
    }
    else {
      wave.boundary = 1
    }
    wave.refreshConstants()
  })
pullDownMenuBoundaryCondition.add("Constant")
pullDownMenuBoundaryCondition.add("Free")

var mousedown = false
cv.element.addEventListener("mousedown", onMouseDownCanvas, false)
cv.element.addEventListener("mousemove", onMouseMoveCanvas, false)
cv.element.addEventListener("mouseup", onMouseUpCanvas, false)
cv.element.addEventListener("mouseout", onMouseOutCanvas, false)

animate()
