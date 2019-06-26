import json
import matplotlib.pyplot as pyplot
import numpy
import scipy.signal as signal

def plot_groupdelay(prefix, irs, delta_min, omega_max):
    pyplot.figure(figsize=(8, 4.8), dpi=100)
    pyplot.title(f"{prefix} Group Delay, N={len(irs[0])}")

    ax = pyplot.gca()
    ax.grid(which="both")

    cmap = pyplot.get_cmap("plasma")
    for index, ir in reversed(list(enumerate(irs))):
        w, gd = signal.group_delay((ir, [1]))
        ax.plot(
            w,
            gd,
            label=str(delta_min + index / (len(irs) - 1)),
            color=cmap(index / len(irs)),
        )

    ax.axvline(
        omega_max * numpy.pi,
        label="$\omega_{\max}$",
        color="black",
        alpha=0.25,
        ls="--",
    )

    ax.set_ylabel("Group delay [samples]")
    ax.set_xlabel("Frequency [rad/sample]")
    ax.legend(
        loc="center right",
        bbox_to_anchor=(1.12, 0.5),
        ncol=1,
        fancybox=True,
        shadow=True,
    )

    offset = 0.3
    ax.set_ylim(delta_min - offset, delta_min + 1 + offset)
    pyplot.savefig(f"img/{prefix}_groupdelay.png", dpi=100)
    # pyplot.plot()

def plot_groupdelay_error(prefix, irs, delta_min, omega_max):
    omega_max *= numpy.pi

    pyplot.figure(figsize=(8, 4.8), dpi=100)
    pyplot.title(f"{prefix} Group Delay, N={len(irs[0])}")

    ax = pyplot.gca()
    ax.grid()

    cmap = pyplot.get_cmap("plasma")

    for index, ir in enumerate(irs):
        w, gd = signal.group_delay((ir, [1]))
        trim = numpy.searchsorted(w, omega_max)
        ax.plot(
            w[:trim + 1],
            gd[:trim + 1] - delta_min - index / (len(irs) - 1),
            label=str(delta_min + index / (len(irs) - 1)),
            color=cmap(index / len(irs)),
            lw=1,
        )

    ax.axvline(omega_max, label="$\omega_{\max}$", color="black", alpha=0.25, ls="--")

    ax.set_ylabel("Error [samples]")
    ax.set_xlabel("Frequency [rad/sample]")
    ax.legend(ncol=2)
    pyplot.savefig(f"img/{prefix}_error.png", dpi=100)
    # pyplot.plot()

def calc_error(irs, delta_min, omega_max):
    """
    irs[fraction][tap]
    """
    error = []
    for index, ir in enumerate(irs):
        w, gd = signal.group_delay((ir, [1]))
        trim = numpy.searchsorted(w, omega_max * numpy.pi)
        error.append(gd[:trim + 1] - delta_min - index / (len(irs) - 1))
    return numpy.sum(numpy.abs(numpy.array(error))) / omega_max / len(irs)

def load_lagrange_ir(order):
    with open("lagrange_ir.json", "r") as infile:
        return json.load(infile)[f"order{order:02d}"]

def load_socp_ir(n_taps):
    with open("table_n_taps.json", "r") as infile:
        table_data = json.load(infile)
    key = str(n_taps)
    return (
        numpy.array(table_data[key]["table"]),
        table_data[key]["delta_min"],
        table_data[key]["omega_max"],
    )

def plot(prefix, irs, delta_min, omega_max):
    print(prefix, calc_error(irs, delta_min, omega_max))
    plot_groupdelay(prefix, irs, delta_min, omega_max)
    # plot_groupdelay_error(prefix, irs, delta_min, omega_max)

order = 32

irs, delta_min, omega_max = load_socp_ir(order)
plot("SOCP", irs, delta_min, omega_max)

irs = load_lagrange_ir(order - 1)
delta_min = int(len(irs[0]) / 2) - 1
plot("Lagrange", irs, delta_min, omega_max)
