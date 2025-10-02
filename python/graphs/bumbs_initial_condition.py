import numpy as np
import matplotlib.pyplot as plt

def my_function(y_var):
    if y_var < 0.5:
        return 25.6 * y_var**4 - 25.6 * y_var**3 + 6.4 * y_var**2 + 0.1
    elif y_var >= 0.5:
        return 25.6 * y_var**4 - 76.8 * y_var**3 + 83.2 * y_var**2 - 38.4 * y_var + 6.5

dim = [0, 1, 0, 1]
param_values = [1]
x_name = r"$x$"
y_name = r"$y$"

plt.rcParams['text.usetex'] = True

fig, ax = plt.subplots(figsize=(8, 6))

y = np.linspace(dim[2], dim[3], 1000)
x = np.zeros(y.size)
for i, y_var in enumerate(y):
    x[i]= my_function(y_var)
ax.plot(x, y, color='black', linewidth=1.5)

ax.text(0.02, 0.5, r"$\Omega_\alpha$", fontsize=20)
ax.text(0.5, 0.5, r"$\Omega_\beta$", fontsize=20)

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
y_ticks = [0, 1]
x_ticks = [0, 1]
ax.set_yticks(y_ticks)
ax.set_xticks(x_ticks)
ax.tick_params(colors='black')

# Bez mřížky
ax.grid(False)

ax.set_aspect("equal")

plt.tight_layout()
plt.savefig("fig\\bumbs_init.png", dpi=300)
