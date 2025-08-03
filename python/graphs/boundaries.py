import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.rcParams['text.usetex'] = True

# Parametry obdélníka (posunutý od os)
x0, y0 = 1.5, 1.5
width, height = 5, 3

# Vytvoření figure a axes
fig, ax = plt.subplots(figsize=(8, 5))

# Přidání obdélníka
rect = patches.Rectangle((x0, y0), width, height, linewidth=2,
                         edgecolor='black', facecolor='none')
ax.add_patch(rect)

# Popisky ke stranám obdélníku
ax.text(x0 + width / 2, y0 - 0.2, r"$\partial \Omega_x^1$", ha='center', va='top', fontsize=22)
ax.text(x0 + width / 2, y0 + height + 0.2, r"$\partial \Omega_x^2$", ha='center', va='bottom', fontsize=22)
ax.text(x0 - 0.2, y0 + height / 2, r"$\partial \Omega_y^1$", ha='right', va='center', fontsize=22)
ax.text(x0 + width + 0.2, y0 + height / 2, r"$\partial \Omega_y^2$", ha='left', va='center', fontsize=22)
ax.text(x0 + width/2, y0 + height / 2, r"$\Omega$", ha='center', va='center', fontsize=22)


# Rozsahy os
ax.set_xlim(0, x0 + width + 1)
ax.set_ylim(0, y0 + height + 1)
ax.set_aspect('equal')

# Skrytí klasických os
for spine in ['left', 'bottom', 'right', 'top']:
    ax.spines[spine].set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

# Vykreslení os jako šipek
arrow_kwargs = dict(arrowstyle="->", color='black', linewidth=1.5)
ax.annotate('', xy=(x0 + width + 0.8, 0), xytext=(-0.5, 0), arrowprops=arrow_kwargs)  # osa x
ax.annotate('', xy=(0, y0 + height + 0.8), xytext=(0, -0.5), arrowprops=arrow_kwargs)  # osa y

# Popisky os
ax.text(x0 + width + 0.9, -0.1, r"$x$", fontsize=22, ha='left', va='center')
ax.text(-0.2, y0 + height + 0.9, r"$y$", fontsize=22, ha='center', va='bottom')

# Uložení
plt.tight_layout()
plt.savefig("fig\\Obdélníková oblast.png", dpi=300)
