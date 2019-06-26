"""
# 使用例

```bash
python3 socp.py -p 31 11 _ 0.9 4
```

# 説明
Putnam and Smith, "Design of Fractional Delay Filters Using Convex Optimization", 1997

cvxopt の matrix は Python のリストと NumPy の ndarray が入力されたときで
行と列の順番が異なる。ドキュメントの例で matrix([[0, 1], [2, 3]]) などとあるときは注意。

```
>>> import numpy
>>> from cvxopt import matrix, solvers
>>>
>>> a = numpy.arange(6).reshape(2, 3)
>>> print(matrix(a))
[ 0  1  2]
[ 3  4  5]

>>> print(matrix(a.tolist()))
[ 0  3]
[ 1  4]
[ 2  5]

>>> print(a)
[[0 1 2]
 [3 4 5]]
>>> print(a.tolist())
[[0, 1, 2], [3, 4, 5]]

```
"""

import argparse
import json
import numpy

from cvxopt import matrix, solvers

def get_b(omega, fraction):
    H_d = numpy.exp(-1j * fraction * omega)
    return numpy.vstack((H_d.real, H_d.imag)).T

def get_A_tilde(n_taps, omega):
    index = numpy.arange(n_taps)
    a = numpy.exp(-1j * index.reshape(1, -1) * omega.reshape(-1, 1))

    A_tilde = []
    for a_i in a:
        a_i0 = numpy.append(a_i, 0)
        A_tilde.append([a_i0.real, a_i0.imag])
    return numpy.array(A_tilde)

def solve(A, b):
    c = numpy.zeros(A.shape[2])
    c[-1] = 1

    rhs = numpy.zeros_like(A[0][0])
    rhs[-1] = 1

    G = [matrix(-numpy.vstack((rhs.reshape(1, -1), A_i))) for A_i in A]
    h = [matrix(numpy.append(0, b_i)) for b_i in b]

    solvers.options["show_progress"] = False
    sol = solvers.socp(matrix(c), Gq=G, hq=h)
    return numpy.array(sol["x"]).flatten()

def createTable(n_taps, n_fraction, delta_min=None, omega_max=0.9, omega_density=4):
    if n_taps < 2:
        raise ValueError("n_taps must be greater than or equal to 2.")
    if n_fraction < 1:
        raise ValueError("n_fraction must be greater than or equal to 1.")
    if delta_min is not None and (delta_min < 0 or delta_min >= n_taps):
        raise ValueError("delta_min must be in range of [0, n_taps - 1].")
    if omega_max < 0 or omega_max > 1:
        raise ValueError("omega_max must be in range of [0, 1].")
    if omega_density < 1:
        raise ValueError("omega_density must be greater than or equal to 1.")

    omega = numpy.linspace(
        0,
        numpy.pi * omega_max,
        int(omega_density * n_taps),
    )

    if delta_min is None:
        delta_min = numpy.floor((n_taps - 1) / 2)

    # fraction is Δ in paper.
    table = [
        -solve(get_A_tilde(n_taps, omega), get_b(omega, fraction))[0:-1]
        for fraction in numpy.linspace(delta_min, delta_min + 1, n_fraction)
    ]
    return {
        "delta_min": delta_min,
        "omega_max": omega_max,
        "table": numpy.array(table),
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--param",
        nargs="+",
        help="[n_tap, n_fraction, delta_min, omega_max, omega_density]",
        default=[32, 11, "_", 0.9, 4],
    )
    args = parser.parse_args()

    table = createTable(
        int(args.param[0]),
        int(args.param[1]),
        int(args.param[2]) if args.param[2].isdigit() else None,
        float(args.param[3]),
        float(args.param[4]),
    )

    table["table"] = table["table"].tolist()
    with open("table.json", "w") as outfile:
        json.dump(table, outfile)
