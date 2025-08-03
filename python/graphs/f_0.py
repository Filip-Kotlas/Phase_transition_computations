import numpy as np
import matplotlib.pyplot as plt

def my_function(x, param):
    return param * x* (1-x)*(x-0.5)

dim = [-0.1, 1.1, -0.06, 0.06]
param_values = [1]
yticks = [-0.05, -0.025, 0, 0.025, 0.05]
x_name = r"$\phi$"
y_name = r"$f_0(\phi)$"

plt.rcParams['text.usetex'] = True

x = np.linspace(dim[0], dim[1], 1000)
fig, ax = plt.subplots(figsize=(8, 6))

for param in param_values:
    y = my_function(x, param)
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

ytick_labels = [f"${y}a$" if y != 0 else "$0$" for y in yticks]
ax.set_yticks(yticks)
ax.set_yticklabels(ytick_labels)
ax.tick_params(labelsize=17)


# Černé hodnoty na osách
ax.tick_params(colors='black')

# Bez mřížky
ax.grid(True)

ax.set_aspect("auto")

plt.tight_layout()
plt.savefig("fig\\f_0.png", dpi=300)
