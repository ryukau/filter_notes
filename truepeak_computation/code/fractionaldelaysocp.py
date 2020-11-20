import json
import numpy

from cvxopt import matrix, solvers

def get_b(omega, delta):
    H_d = numpy.exp(-1j * delta * omega)
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

def createTable(n_taps, n_fraction, delta_min=None, omega_max=0.9, omega_density=1):
    omega = numpy.linspace(
        0,
        numpy.pi * omega_max,
        omega_density * n_taps,
    )

    if delta_min is None:
        delta_min = n_taps / 2 - 1

    table = [
        -solve(get_A_tilde(n_taps, omega), get_b(omega, delta))[0:-1]
        for delta in numpy.linspace(delta_min, delta_min + 1, n_fraction)
    ]
    return {
        "delta_min": delta_min,
        "omega_max": omega_max,
        "table": numpy.array(table),
    }

def printCpp():
    with open("socp.json", "r", encoding="utf-8") as fi:
        data = json.load(fi)

    table = data["table"]
    if len(table[0]) % 2 == 0:  # length of the filter is even.
        table = table[1:-1]
    else:
        table.pop((len(table) - 1) // 2)
        table = table[:-1]

    for tbl in table:
        text = "{"
        for idx, value in enumerate(tbl):
            if idx != len(tbl) - 1:
                text += f"Sample({value}), "
            else:
                text += f"Sample({value})"
        print(text + "},")

if __name__ == "__main__":
    table = createTable(7, 5, omega_max=0.65, omega_density=1)

    table["table"] = table["table"].tolist()
    with open("socp.json", "w") as outfile:
        json.dump(table, outfile)

    print()
    printCpp()
