import json
import numpy
import scipy.signal as signal
import matplotlib.ticker
import matplotlib.pyplot as pyplot

def get_filter_coefficients(table_data, delta):
    if delta < 0 or delta > 1:
        return None
    elif delta == 1:
        return numpy.array(table_data["table"])[-1]
    table = numpy.array(table_data["table"])
    length = table.shape[0] - 1
    pos = delta * length
    low = int(pos)
    return table[low] + (pos - low) * (table[low + 1] - table[low])

with open("table.json", "r") as infile:
    table_data = json.load(infile)

source = signal.unit_impulse(33)

pyplot.figure(figsize=(8, 4.8), dpi=100)
pyplot.title(f"Impulse Response, N={len(table_data['table'][0])}")
cmap = pyplot.get_cmap("plasma")
ax = pyplot.gca()
ax.plot(source, label="source", color="black")
for fraction in numpy.linspace(0, 1, 11):
    fir = get_filter_coefficients(table_data, fraction)
    dest = signal.lfilter(fir, [1], source)
    ax.plot(
        dest,
        label=f"Δ={table_data['delta_min'] + fraction}",
        alpha=0.8,
        color=cmap(fraction * 0.9999),
    )
ax.grid(which="both")
ax.legend()
ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(base=1.0))
# pyplot.savefig("img/impulse_response.png", dpi=100)
pyplot.show()
