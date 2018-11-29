# Energy Decay Curveの実装

$$
EDC(t) \triangleq \int_{t}^{\infty} h^2(\tau) d\tau
$$

$\triangleq$ は定義するという意味。

$h(\tau)$ はシステムのインパルス応答。

```javascript
function getEDC(ir) {
  var edc = new Array(ir.length)
  var last = ir.length - 1
  edc[last] = ir[last] * ir[last]
  for (var t = ir.length - 2; t >= 0; --t) {
    edc[t] = ir[t] * ir[t] + edc[t + 1]
  }
  return edc
}

var ir = generateIR() // インパルス応答。
var edc = getEDC(ir)
```

```python
import numpy
from scipy.signal import *

def generateIR():
    pass

ir = generateIR()
edc = numpy.cumsum((ir * ir)[::-1])[::-1]
```
