import numpy as np
import matplotlib.pyplot as plt
import json


with open("phase.json", "r", encoding="utf-8") as fi:
    data = json.load(fi)

cmap = plt.get_cmap("viridis")
for index, (key, value) in enumerate(data.items()):
    plt.plot(value, color=cmap(index / len(data)), label=key)
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
