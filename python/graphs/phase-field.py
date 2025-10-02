import numpy as np
import matplotlib.pyplot as plt

def my_function(x_var):
    return 0.5 * np.tanh(-3/0.1 * (x_var - 0.5)) + 0.5

x_name = r"$x$"
y_name = r"$\phi$"

plt.rcParams['text.usetex'] = True

x = np.linspace(-0.1, 1.3, 1000)
y = my_function(x)

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(x, y, color='black', linewidth=1.5)

ax.axvline(x=0.4, color='gray', linestyle='--', linewidth=1)
ax.axvline(x=0.6, color='gray', linestyle='--', linewidth=1)
ax.axhline(y=0.5, color='gray', linestyle='--', linewidth=1)

# Dvoustranná šipka mezi x=0.4 a x=0.6 nad křivkou
arrow_y = 1.05  # výška, kam šipku umístit – uprav dle potřeby
ax.annotate(
    "", 
    xy=(0.6, arrow_y), 
    xytext=(0.4, arrow_y), 
    arrowprops=dict(arrowstyle='<->', color='black', linewidth=1.5)
)

# Popisky oblastí nad grafem
label_y = arrow_y + 0.16  # výška nad šipkou – uprav dle potřeby

ax.text(0.2, label_y, r"$\Omega_\alpha$", ha='center', va='bottom', fontsize=18)
ax.text(0.5, label_y, r"$\Omega_\Gamma$", ha='center', va='bottom', fontsize=18)
ax.text(0.8, label_y, r"$\Omega_\beta$", ha='center', va='bottom', fontsize=18)

# Popisek nad šipku (volitelné)
ax.text(0.5, arrow_y + 0.025, r"$\xi$", ha='center', va='bottom', fontsize=18)

# Rozsah os
ax.set_xlim(-0.1, 1.1)
ax.set_ylim(-0.1, 1.2)

# Zapnutí klasických os a jejich vzhledu
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_color('black')

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax.tick_params(labelsize=17, colors='black')

# Popisky os
ax.text(1.02, -0.02, x_name, ha='right', va='top', fontsize=20, transform=ax.transAxes)
ax.text(-0.03, 1.02, y_name, ha='left', va='bottom', fontsize=20, transform=ax.transAxes)

ax.set_xticks([0, 0.5, 1.0])
ax.set_yticks([0, 0.5, 1])

# Poměr osy
ax.set_aspect('auto')

plt.tight_layout()
plt.savefig("fig\\phase-field.png", dpi=300)
