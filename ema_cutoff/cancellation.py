import numpy as np
from mpmath import mp


def example_pi():
    dtype = np.float32
    mp.dps = 50

    x_mp = mp.mpf(355) / mp.mpf(113)
    y_mp = dtype(mp.pi - x_mp)

    x_np = dtype(355.0) / dtype(113.0)
    y_np = dtype(np.pi) - x_np

    rel = np.abs((y_mp - y_np) / y_mp, dtype=dtype)
    ulp = np.abs(y_mp - y_np, dtype=dtype) / np.spacing(y_mp, dtype=dtype)

    text = f"""
π (mp)  : {mp.pi:.9e}
π (f32) : {dtype(np.pi):.9e}

355 / 113 (mp)  : {x_mp:.9e}
355 / 113 (f32) : {x_np:.9e}

π - 355 / 113 (mp)  : {y_mp:.9e}
π - 355 / 113 (f32) : {y_np:.9e}

Rel Error : {rel: .2e}
ULP Error : {ulp: .2e}
"""
    print(text)


def example_1mcos():
    dtype = np.float32
    mp.dps = 50

    x = dtype(1e-1)

    y_mp = dtype(mp.mpf(1) - mp.cos(mp.mpf(float(x))))
    y_np = dtype(1 - np.cos(x))

    rel = np.abs((y_mp - y_np) / y_mp, dtype=dtype)
    ulp = np.abs(y_mp - y_np, dtype=dtype) / np.spacing(y_mp, dtype=dtype)

    text = f"""
x : {x:.9e}

1 - cos(x) (mp)  : {y_mp:.9e}
1 - cos(x) (f32) : {y_np:.9e}

Rel Error : {rel: .2e}
ULP Error : {ulp: .2e}
"""
    print(text)


example_pi()
# example_1mcos()
