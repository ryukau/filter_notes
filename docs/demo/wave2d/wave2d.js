/* math functions */
class Solver {
  // Ax = b を解く。
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

    this.aDiagonal = new Array(a.length)
    for (var i = 0; i < a.length; ++i) {
      this.aDiagonal[i] = a[i][i]
    }

    this.x_prev = new Array(a.length).fill(0)
  }

  setA(aReduced, aDiagonal) {
    this.aReduced = aReduced
    this.aDiagonal = aDiagonal
    this.x_prev = new Array(aDiagonal.length).fill(0)
  }

  isInTolerance(x) {
    for (var i = 0; i < x.length; ++i) {
      if (Math.abs(x.array[i].value - this.x_prev[i]) > this.tolerance) {
        return false
      }
    }
    return true
  }

  solveBXNode(b, x) {
    for (var i = 0; i < x.length; ++i) {
      x.array[i].value = 0
    }

    for (var iter = 0; iter < this.maxIteration; ++iter) {
      for (var i = 0; i < x.length; ++i) {
        this.x_prev[i] = x.array[i].value
      }

      for (var i = 0; i < x.length; ++i) {
        var sum = b.array[i].value
        for (var j = 0; j < this.aReduced[i].length; ++j) {
          sum -= this.aReduced[i][j][0] * x.array[this.aReduced[i][j][1]].value
        }
        x.array[i].value = sum / this.aDiagonal[i]
      }

      for (var i = 0; i < x.length; ++i) {
        this.x_prev[i] = x.array[i].value
      }

      for (var i = x.length - 1; i >= 0; --i) {
        var sum = b.array[i].value
        for (var j = 0; j < this.aReduced[i].length; ++j) {
          sum -= this.aReduced[i][j][0] * x.array[this.aReduced[i][j][1]].value
        }
        x.array[i].value = sum / this.aDiagonal[i]
      }

      if (this.isInTolerance(x)) break
    }
  }

  solveBXNodeJacobi(b, x) {
    for (var i = 0; i < x.length; ++i) {
      x.array[i].value = 0
    }

    var omega = 2 / 3
    var omega1 = 1 - omega

    for (var iter = 0; iter < this.maxIteration; ++iter) {
      for (var i = 0; i < x.length; ++i) {
        this.x_prev[i] = x.array[i].value
      }

      for (var i = 0; i < x.length; ++i) {
        var sum = b.array[i].value
        for (var j = 0; j < this.aReduced[i].length; ++j) {
          sum -= this.aReduced[i][j][0] * this.x_prev[this.aReduced[i][j][1]]
        }
        x.array[i].value
          = omega * sum / this.aDiagonal[i]
          + omega1 * this.x_prev[i]
      }

      for (var i = 0; i < x.length; ++i) {
        this.x_prev[i] = x.array[i].value
      }

      for (var i = x.length - 1; i >= 0; --i) {
        var sum = b.array[i].value
        for (var j = 0; j < this.aReduced[i].length; ++j) {
          sum -= this.aReduced[i][j][0] * this.x_prev[this.aReduced[i][j][1]]
        }
        x.array[i].value
          = omega * sum / this.aDiagonal[i]
          + omega1 * this.x_prev[i]
      }

      if (this.isInTolerance(x)) break
    }
  }
}

/* wave */
class Array2D {
  constructor(width, height) {
    this.width = width
    this.height = height
    this.length = width * height

    this.lattice = new Array(this.width)
    this.array = new Array(this.length)
    for (var x = 0; x < this.width; ++x) {
      this.lattice[x] = new Array(this.height)
      for (var y = 0; y < this.height; ++y) {
        var node = { value: 0 }
        this.lattice[x][y] = node
        this.array[x + y * this.width] = node
      }
    }
  }
}

class DampedWave2DNewmarkBeta {
  constructor(width, height, c, dx, dt, a, k) {
    this.width = width
    this.height = height
    this.length = width * height
    this.xLast = width - 1
    this.yLast = height - 1
    this.initArrays()

    this.c = c
    this.dt = dt
    this.dx = dx
    this.a = a
    this.k = k

    this.beta = 1 / 4

    this.boundary = 1

    this.solver = new Solver()
    this.solver.tolerance = 1e-9
    this.solver.maxIteration = 8
    this.vectorB = new Array2D(this.width, this.height)

    this.refreshConstants()

    this.isPicked = false
    this.pickX = 0
    this.pickY = 0
  }

  initArrays() {
    this.wave = new Array2D(this.width, this.height)
    this.velocity = new Array2D(this.width, this.height)

    this.acceleration = []
    for (var i = 0; i < 2; ++i) {
      this.acceleration.push(new Array2D(this.width, this.height))
    }
  }

  refreshConstants() {
    var dt2 = this.dt * this.dt

    this.C2 = this.c * this.c / dt2
    this.C3 = dt2 * (0.5 - this.beta)
    this.C7 = this.dt / 2
    this.C8 = dt2 * this.beta
    this.C1 = - this.C2 * this.C8
    this.C6 = this.k + 4 * this.C2
    this.C0 = 1 + this.a * this.C7 + this.C6 * this.C8
    this.C4 = this.C3 * this.C6 + this.a * this.C7
    this.C5 = this.a + this.dt * this.C6

    this.initMatrix()
  }

  initMatrix() {
    // aReduced = [[value, index], ...]
    var aReduced = new Array(this.length)
    for (var i = 0; i < aReduced.length; ++i) {
      aReduced[i] = []
    }

    for (var x = 0; x < this.width; ++x) {
      for (var y = 0; y < this.height; ++y) {
        var index = x + y * this.width
        if (x === 0) {
          aReduced[index].push([this.boundary * this.C1, index + 1])
        }
        else if (x === this.xLast) {
          aReduced[index].push([this.boundary * this.C1, index - 1])
        }
        else {
          aReduced[index].push([this.C1, index - 1])
          aReduced[index].push([this.C1, index + 1])
        }
        if (y === 0) {
          aReduced[index].push([this.boundary * this.C1, index + this.width])
        }
        else if (y === this.yLast) {
          aReduced[index].push([this.boundary * this.C1, index - this.width])
        }
        else {
          aReduced[index].push([this.C1, index - this.width])
          aReduced[index].push([this.C1, index + this.width])
        }
      }
    }

    var aDiagonal = new Array(this.length).fill(this.C0)

    this.solver.setA(aReduced, aDiagonal)
  }

  sumNeighbor(array, x, y) {
    const sumX = x === 0
      ? this.boundary * array.lattice[x + 1][y].value
      : x === this.xLast
        ? this.boundary * array.lattice[x - 1][y].value
        : array.lattice[x - 1][y].value + array.lattice[x + 1][y].value
    const sumY = y === 0
      ? this.boundary * array.lattice[x][y + 1].value
      : y === this.yLast
        ? this.boundary * array.lattice[x][y - 1].value
        : array.lattice[x][y - 1].value + array.lattice[x][y + 1].value
    return sumX + sumY
  }

  step() {
    this.acceleration.unshift(this.acceleration.pop())

    var last = this.length - 1

    for (var y = 0; y < this.width; ++y) {
      for (var x = 0; x < this.height; ++x) {
        this.vectorB.lattice[x][y].value
          = this.C2 * (
            this.sumNeighbor(this.wave, x, y)
            + this.dt * this.sumNeighbor(this.velocity, x, y)
            + this.C3 * this.sumNeighbor(this.acceleration[1], x, y)
          )
          - this.C4 * this.acceleration[1].lattice[x][y].value
          - this.C5 * this.velocity.lattice[x][y].value
          - this.C6 * this.wave.lattice[x][y].value
      }
    }

    this.solver.solveBXNode(this.vectorB, this.acceleration[0])
    // this.solver.solveBXNodeJacobi(this.vectorB, this.acceleration[0])

    for (var i = 0; i < this.length; ++i) {
      this.wave.array[i].value += this.dt * this.velocity.array[i].value
        + this.C3 * this.acceleration[1].array[i].value
        + this.C8 * this.acceleration[0].array[i].value
      this.velocity.array[i].value += this.C7 * (
        this.acceleration[1].array[i].value
        + this.acceleration[0].array[i].value
      )
    }

    if (this.isPicked) {
      this.acceleration[0].lattice[this.pickX][this.pickY].value = 4 * this.c
    }
  }

  reset() {
    for (var i = 0; i < this.length; ++i) {
      this.wave.array[i].value = 0
      this.velocity.array[i].value = 0
      this.acceleration[0].array[i].value = 0
      this.acceleration[1].array[i].value = 0
    }
  }

  pick(x, y) {
    this.isPicked = x >= 0 && x < this.width && y >= 0 && y < this.height
    if (this.isPicked) {
      this.pickX = x
      this.pickY = y
    }
  }
}

/* UI */
function pickWave(event) {
  var rect = event.target.getBoundingClientRect()
  var x = event.clientX - rect.left
  var y = event.clientY - rect.top

  wave.pick(
    Math.floor(waveSize * x / cv.width),
    Math.floor(waveSize * y / cv.height)
  )
}

function onMouseDownCanvas(event) {
  mousedown = true
  pickWave(event)
}

function onMouseMoveCanvas(event) {
  if (!mousedown) {
    return
  }
  pickWave(event)
}

function onMouseUpCanvas(event) {
  mousedown = false
  wave.pick(-1, -1)
}

function onMouseOutCanvas(event) {
  onMouseUpCanvas(event)
}

function updateCanvas() {
  cv.clearWhite()

  var xCenter = cv.width / 2
  var yCenter = cv.height / 2

  var position = wave.wave
  var pixels = cv.pixels
  for (var x = 0; x < cv.width; ++x) {
    for (var y = 0; y < cv.height; ++y) {
      var index = (y * width + x) * 4
      var ix = Math.floor(waveSize * x / cv.width)
      var iy = Math.floor(waveSize * y / cv.width)
      var color = Math.floor((1 + position.lattice[ix][iy].value) * 127)
      if (color < 0) {
        pixels[index + 0] = 0
        pixels[index + 1] = 255
        pixels[index + 2] = 0
        pixels[index + 3] = 255
      }
      else if (color >= 256) {
        pixels[index + 0] = 255
        pixels[index + 1] = 0
        pixels[index + 2] = 0
        pixels[index + 3] = 255
      }
      else {
        pixels[index + 0] = color
        pixels[index + 1] = color
        pixels[index + 2] = color
        pixels[index + 3] = 255
      }
    }
  }
  cv.putPixels()
}

function animate() {
  updateCanvas()
  wave.step()
  if (isAnimating) {
    requestAnimationFrame(animate)
  }
}

/* main */
var width = 512
var height = 512
var waveSize = 128

var isAnimating = true

var wave = new DampedWave2DNewmarkBeta(
  waveSize, waveSize, 4, 0.1, 100 / 60, 1, 1)

var cv = new Canvas(document.body, width, height)
var divWave = new Div(document.body, "divWave")
var buttonReset = new Button(divWave.element,
  "Reset",
  () => { wave.reset() })
var buttonStep = new Button(divWave.element,
  "pause",
  () => {
    isAnimating = !isAnimating
    buttonStep.element.value = isAnimating ? "pause" : "animate"
    animate()
  })
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
