import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
#from scipy.spatial import cKDTree

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

# Funkce pro nalezení bodů na hranici
def find_boundary(x, y, values, threshold=0.5):
    # Předpokládáme, že body jsou na mřížce (nebo téměř na mřížce)
    grid_shape = (len(np.unique(x)), len(np.unique(y)))
    
    # Převod na matici (2D mřížku)
    grid_x = x.reshape(grid_shape)
    grid_y = y.reshape(grid_shape)
    grid_values = values.reshape(grid_shape)
    
    # Vyhledání bodů na hranici
    boundary_x, boundary_y = [], []
    offsets = [-1, 0, 1, 0]
    for i in range(grid_values.shape[0] - 1):
        for j in range(grid_values.shape[1] - 1):
            # Kontrola hranic v obou osách
            neighbour_over = False
            neighbour_below = False
            for k in range(4):
                if grid_values[i + offsets[k], j + offsets[(k+1)%4]] >= threshold:
                    neighbour_over = True
                if grid_values[i + offsets[k], j + offsets[(k+1)%4]] <= threshold:
                    neighbour_below = True

            if neighbour_over and neighbour_below and grid_values[i, j] >= threshold:
                boundary_x.append(grid_x[i, j])
                boundary_y.append(grid_y[i, j])
    
    return np.array(boundary_x), np.array(boundary_y)

def get_analytic_solution(x, y, r0, t, eps):
    if r0**2 + 2*t < 0:
        return np.array([]), np.array([])  # Žádný kruh neexistuje
    r = np.sqrt(r0**2 + 2*t)

    # Vyhledání bodů na hranici
    boundary_x, boundary_y = [], []

    for xi, yi in zip(x, y):
        if abs(np.sqrt(xi**2 + yi**2) - r) < eps:
            boundary_x.append(xi)
            boundary_y.append(yi)
    #return np.array([]), np.array([])
    return np.array(boundary_x), np.array(boundary_y)

# Inicializace grafu
fig, ax = plt.subplots()
sc = ax.scatter([], [], c=[], cmap='viridis', s=1)
boundary_sc = ax.scatter([], [], color='red', s=1, label='Numerical solution')
analytic_solution_sc = ax.scatter([], [], color='green', s=1, label='Analytic solution')
ax.set_xlim(np.loadtxt(os.path.join(folder_path, files[0]))[:, 0].min(), 
            np.loadtxt(os.path.join(folder_path, files[0]))[:, 0].max())
ax.set_ylim(np.loadtxt(os.path.join(folder_path, files[0]))[:, 1].min(), 
            np.loadtxt(os.path.join(folder_path, files[0]))[:, 1].max())
ax.legend()
ax.set_aspect('equal')
title = ax.set_title("Frame: 0")

# Funkce pro aktualizaci snímku v animaci
def update(frame):
    x, y, values = load_data(os.path.join(folder_path, files[frame]))
    sc.set_offsets(np.column_stack((x, y)))
    sc.set_array(values)

    boundary_x, boundary_y = find_boundary(x, y, values, threshold=0.5)
    boundary_sc.set_offsets(np.column_stack((boundary_x, boundary_y)))
    analytic_boundary_x, analytic_boundary_y = get_analytic_solution(x, y, 0.5, frame*0.001, (x[1]-x[0])/2)
    analytic_solution_sc.set_offsets(np.column_stack((analytic_boundary_x, analytic_boundary_y)))

    title.set_text(f"Frame: {frame}")
    print("Updating frame: ", frame)
    #return boundary_sc, analytic_solution_sc, title #pridat sc kdybych chtěl vykreslovat celou tu funkci
    return boundary_sc, title #pridat sc kdybych chtěl vykreslovat celou tu funkci
    #return sc

# Vytvoření animace
ani = animation.FuncAnimation(fig, update, frames=len(files), blit=True)

# Zobrazení animace
#plt.colorbar(sc)

output_file = 'ACE.gif'
ani.save(output_file, writer='pillow', fps=10)
print("Finished. Saved as:", output_file)