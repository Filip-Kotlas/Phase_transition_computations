import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Nastavte cestu k adresáři, kde jsou soubory
folder_path = 'Results'

# Získání seznamu všech souborů
files = sorted([f for f in os.listdir(folder_path) if f.endswith('.txt')])

# Předpokládám, že všechny soubory mají stejný formát a stejný rozsah
def load_data(file_path):
    # Načtení dat z jednoho souboru
    data = np.loadtxt(file_path)
    x = data[:, 0]
    y = data[:, 1]
    values = data[:, 2]
    return x, y, values

# Inicializace grafu
fig, ax = plt.subplots()
sc = ax.scatter([], [], c=[], cmap='viridis', s=5)
"""
ax.set_xlim(np.min([np.loadtxt(os.path.join(folder_path, f))[:, 0].min() for f in files]), 
            np.max([np.loadtxt(os.path.join(folder_path, f))[:, 0].max() for f in files]))
ax.set_ylim(np.min([np.loadtxt(os.path.join(folder_path, f))[:, 1].min() for f in files]), 
            np.max([np.loadtxt(os.path.join(folder_path, f))[:, 1].max() for f in files]))
"""
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.set_aspect('equal')
title = ax.set_title("Frame: 0")

# Funkce pro aktualizaci snímku v animaci
def update(frame):
    x, y, values = load_data(os.path.join(folder_path, files[frame]))
    sc.set_offsets(np.column_stack((x, y)))
    sc.set_array(values)
    title.set_text(f"Frame: {frame}")
    print("Updating frame: ", frame)
    return sc,

# Vytvoření animace
ani = animation.FuncAnimation(fig, update, frames=len(files), interval=1, blit=True)

# Zobrazení animace
plt.colorbar(sc)

ani.save('ACE.gif', writer='pillow', fps=10)
print("Finished")