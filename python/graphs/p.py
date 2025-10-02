import numpy as np
import matplotlib.pyplot as plt

def my_function(x):
    return x**3*(6*x**2 -15*x + 10)

dim = [-0.15, 1.15, -0.1, 1.1]
param_values = [1]
x_name = r"$\phi$"
y_name = r"$p(\phi)$"

plt.rcParams['text.usetex'] = True

x = np.linspace(dim[0], dim[1], 1000)
fig, ax = plt.subplots(figsize=(8, 6))

for param in param_values:
    y = my_function(x)
    ax.plot(x, y, color='black', linewidth=1.5)


# Nastav rozsah os
ax.set_xlim(dim[0], dim[1])
ax.set_ylim(dim[2], dim[3])

# Zobraz všechny okraje (spines)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')

# Tick pozice klasicky
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')

# Souřadnicové osy
ax.axhline(0, color='black', linewidth=0.5)
ax.axvline(0, color='black', linewidth=0.5)

# Popisky osy – zarovnané k pravému a hornímu okraji
ax.set_xlabel(x_name, color='black', labelpad=10, fontsize=20)
ax.set_ylabel(y_name, color='black', labelpad=20, fontsize=20)
ax.yaxis.label.set_rotation(0)

ax.tick_params(labelsize=17)


# Černé hodnoty na osách
ax.tick_params(colors='black')

# Bez mřížky
ax.grid(True)

ax.set_aspect("auto")

plt.tight_layout()
plt.savefig("fig\\p.png", dpi=300)
