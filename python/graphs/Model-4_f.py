import numpy as np
import matplotlib.pyplot as plt

def my_function(x, param):
    return x*(1-x)*(x-0.5-param*x*(1-x))

dim = [-0.15, 1.15, -0.25, 0.25]
param_values = [-3, -2, -1, 0, 1, 2, 3]
yticks = [-0.20, -0.1,  0, 0.1, 0.2,]
x_name = r"$\phi$"
y_name = r"$f(\phi)$"

cmap = plt.cm.rainbow
norm = plt.Normalize(vmin=min(param_values), vmax=max(param_values))

plt.rcParams['text.usetex'] = True

x = np.linspace(dim[0], dim[1], 1000)
fig, ax = plt.subplots(figsize=(8, 6))

for param in param_values:
    y = my_function(x, param)
    color = cmap(norm(param))
    ax.plot(x, y, color=color, linewidth=1.5, label=fr"$b\xi F(c) = {param}$")

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

ax.legend(loc="upper right", fontsize=13, frameon=True, facecolor='white', edgecolor='black')

ax.set_aspect("auto")

plt.tight_layout()
plt.savefig("fig\\Model-4_f.png", dpi=300)
