/* wave */
class Wave1D {
  constructor(length, c, dx, dt, attenuation) {
    this.length = length
    this.reset()

    this.attenuation = attenuation

    this.c = c
    this.dt = dt
    this.dx = dx

    this.refreshConstants()

    this.edge = 0

    // 0: const
    // 1: free
    this.boundaryCondition = 0
  }

  step() {
    this.wave.unshift(this.wave.pop())

    // 左端はユーザからの入力。
    this.wave[0][0] = wave.edge

    var last = this.wave[0].length - 1

    // 右端。
    if (this.boundaryCondition === 1) {
      // 自由端。
      this.wave[0][last] = this.attenuation * (
        this.alpha * (this.wave[1][last - 1] + this.wave[1][last - 1])
        + this.beta * this.wave[1][last]
        - this.wave[2][last]
      )
    }
    else {
      // 固定端。
      this.wave[0][last] = 0
    }

    for (var x = 1; x < last; ++x) {
      this.wave[0][x] = this.attenuation * (
        this.alpha * (this.wave[1][x + 1] + this.wave[1][x - 1])
        + this.beta * this.wave[1][x]
        - this.wave[2][x]
      )
    }

    // var half = Math.floor(last / 2)
    // this.wave[0][half] = Math.max(Math.min(
    //   this.wave[0][half] + wave.edge / 100, 0.5), -0.5)
  }

  pulse(position, width, height) {
    // Hanning Window without modulo.
    var left = Math.floor(position - width * 0.5)
    var start = index < 0 ? width + left + 1 : 0
    for (var i = start; i < width; ++i) {
      var index = left + i
      if (index >= this.length)
        break
      this.wave[0][index] += height
        * 0.5 * (1 - Math.cos(2 * Math.PI * i / (width - 1)))
    }
  }

  reset() {
    // u(x, t) -> this.wave[t][x]
    this.wave = []
    for (var i = 0; i < 3; ++i) {
      this.wave.push(new Array(this.length).fill(0))
    }
  }

  refreshConstants() {
    this.alpha = (this.c * this.dt / this.dx) ** 2
    this.beta = 2 * (1 - this.alpha)
  }
}

/* UI */
function onMouseDownCanvas(event) {
  mousedown = true
  var rect = event.target.getBoundingClientRect()
  var x = event.clientX - rect.left
  var y = event.clientY - rect.top

  var height = y / cv.height - 0.5
  // wave.pulse(x, width / 4, height / 100)
  wave.edge = height
}

function onMouseMoveCanvas(event) {
  if (!mousedown)
    return
  var rect = event.target.getBoundingClientRect()
  var x = event.clientX - rect.left
  var y = event.clientY - rect.top

  wave.edge = y / cv.height - 0.5
}

function onMouseUpCanvas(event) {
  mousedown = false
  wave.edge = 0
}

function onMouseOutCanvas(event) {
  onMouseUpCanvas(event)
}

function updateCanvas() {
  cv.clearWhite()
  cv.context.save()
  cv.context.translate(0, cv.height / 2)

  cv.context.lineWidth = 10
  cv.context.strokeStyle = "#4488ff"
  cv.context.fillStyle = "#4488ff"
  cv.context.lineJoin = "round"
  cv.context.beginPath()
  cv.context.moveTo(0, wave.wave[0][0] * cv.width / 2)
  for (let x = 1; x < width; ++x) {
    cv.context.lineTo(x, wave.wave[0][x] * cv.width / 2)
  }
  cv.context.lineTo(width, wave.wave[0][width - 1] * cv.width / 2)
  cv.context.lineTo(cv.width, cv.height)
  cv.context.lineTo(0, cv.height)
  cv.context.closePath()
  cv.context.fill()

  var yLine = wave.edge * cv.height
  cv.context.lineWidth = 1
  cv.context.strokeStyle = "#404040"
  cv.context.beginPath()
  cv.context.moveTo(0, yLine)
  cv.context.lineTo(cv.width, yLine)
  cv.context.stroke()

  cv.context.restore()
}

function animate() {
  updateCanvas()
  for (var i = 0; i < 4; ++i)
    wave.step()
  requestAnimationFrame(animate)
}

/* main */
var width = 512
var height = 512

var wave = new Wave1D(width, 4, 0.1, 1 / 60, 1)

var divWave = new Div(document.body, "divWave")
var cv = new Canvas(divWave.element, width, height)
var buttonReset = new Button(divWave.element,
  "Reset",
  () => { wave.reset() })
var numberC = new NumberInput(divWave.element,
  "c", wave.c, 0, 6, 0.1,
  (value) => { wave.c = value; wave.refreshConstants() })
var numberDx = new NumberInput(divWave.element,
  "dx", wave.dx, 0.1, 1, 0.1,
  (value) => { wave.dx = value; wave.refreshConstants() })
var numberAttenuation = new NumberInput(divWave.element,
  "attenuation", wave.attenuation, 0.98, 1, 0.0001,
  (value) => { wave.attenuation = value })
var pullDownMenuBoundaryCondition = new PullDownMenu(divWave.element,
  "Boundary Condition",
  (value) => {
    if (value === "Free") {
      wave.boundaryCondition = 1
    }
    else {
      wave.boundaryCondition = 0
    }
  })
pullDownMenuBoundaryCondition.add("Constant")
pullDownMenuBoundaryCondition.add("Free")
// pullDownMenuBoundaryCondition.add("Ring")

var mousedown = false
cv.element.addEventListener("mousedown", onMouseDownCanvas, false)
cv.element.addEventListener("mousemove", onMouseMoveCanvas, false)
cv.element.addEventListener("mouseup", onMouseUpCanvas, false)
cv.element.addEventListener("mouseout", onMouseOutCanvas, false)

animate()
