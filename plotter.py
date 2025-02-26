from typing import List, Tuple
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


class Plotter():
    """Parent class of specialized plotters."""

    def __init__(self, folder_path: str,
                 draw_function: bool=False,
                 draw_num_boundary: bool=True,
                 draw_anal_boundary: bool=False):
        self.folder = folder_path
        self.file_paths = sorted([os.path.join(self.folder, f) for f in os.listdir(folder_path) if f.endswith('.txt')])
        self.fig = None
        self.draw_function = draw_function
        self.draw_num_boundary = draw_num_boundary
        self.draw_anal_boundary = draw_anal_boundary

    def _load_data(self, file_path: str) -> Tuple[List[float], List[float], List[float]]:
        """
        Loads data from a file.
        
        Args:
            file_path (str): Path of the file with data.

        Returns:
            x: List of "x" coordinates of the data.
            y: List of "y" coordinates of the data.
            values: List of values of the data.
        """
        data = np.loadtxt(file_path)
        x = data[:, 0]
        y = data[:, 1]
        values = data[:, 2]
        return x, y, values

    def reset_draw_settings(self,
                            draw_function: bool=False,
                            draw_num_boundary: bool=True,
                            draw_anal_boundary: bool=False):
        self.draw_function = draw_function
        self.draw_num_boundary = draw_num_boundary
        self.draw_anal_boundary = draw_anal_boundary

    def check_for_visualisation_folder(self):
        path_to_visual = os.path.join(self.folder, "info")
        if not os.path.exists(path_to_visual):
            os.makedirs(path_to_visual)
            print("Creating folder \"info\" in ", self.folder)

    def update(self, frame: int):
        pass

    def save_animation(self) -> None:
        self.check_for_visualisation_folder()
        num_frames = len(self.file_paths)
        ani = animation.FuncAnimation(self.fig, self.update, frames=num_frames, blit=True)
        output_file = os.path.join(self.folder, "info", "ACE.gif")
        ani.save(output_file, writer='pillow', fps=max(num_frames/10, 10))
        print("Finished the animation. Saved at:", output_file)

    def save_frame(self, frame: int) -> None:
        self.check_for_visualisation_folder()
        self.update(frame)
        output_file = os.path.join(self.folder, "info", f"Frame_{frame:05d}.jpg")
        self.fig.savefig(output_file, dpi=300)
        print(f"Saved the frame {frame} at:", output_file)

class ScatterPlotter2D(Plotter):
    """Plotter of 2D scatter plot."""
    def __init__(self,
                 folder_path: str,
                 draw_function: bool=False,
                 draw_num_boundary: bool=True,
                 draw_anal_boundary: bool=False):
        super().__init__(folder_path,
                         draw_function,
                         draw_num_boundary,
                         draw_anal_boundary)

        x, y, _ = self._load_data(self.file_paths[0])
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(x.min(), x.max())
        self.ax.set_ylim(y.min(), y.max())
        self.ax.set_aspect('equal')
        self.title = self.ax.set_title("Frame: 0")

        if self.draw_function:
            self.function_sc = self.ax.scatter([], [], c=[], cmap='viridis', s=1)
            self.function_sc.set_offsets(np.column_stack((x, y)))
        if self.draw_num_boundary:
            self.num_bound_sc = self.ax.scatter([], [], color='red', s=1, label='Num. sol.')
        if self.draw_anal_boundary:
            self.anal_bound_sc = self.ax.scatter([], [], color='green', s=1, label='Anal. sol.')

        self.ax.legend()

    def update(self, frame: int):
        x, y, values = self._load_data(self.file_paths[frame])
        self.title.set_text(f"Frame: {frame}")
        print("Updating frame: ", frame)
        artist = []

        if self.draw_function:
            self.function_sc.set_array(values)
            artist.append(self.function_sc)

        if self.draw_num_boundary:
            num_boundary_x, num_boundary_y = self._find_boundary(x, y, values, threshold=0.5)
            self.num_bound_sc.set_offsets(np.column_stack((num_boundary_x, num_boundary_y)))
            artist.append(self.num_bound_sc)

        if self.draw_anal_boundary:
            analytic_boundary_x, analytic_boundary_y = self.get_analytic_solution(x, y, 0.5, frame*0.001, (x[1]-x[0])/2)
            self.anal_bound_sc.set_offsets(np.column_stack((analytic_boundary_x, analytic_boundary_y)))
            artist.append(self.anal_bound_sc)

        artist.append(self.title)
        return artist

    def _find_boundary(self, x, y, values, threshold=0.5):
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

    def get_analytic_solution(self, x, y, r0, t, eps):
        # There is no circle.
        if r0**2 + 2*t < 0:
            return np.array([]), np.array([])
        r = np.sqrt(r0**2 + 2*t)

        # Locating boundary points
        analytic_boundary_x, analytic_boundary_y = [], []

        for xi, yi in zip(x, y):
            if abs(np.sqrt(xi**2 + yi**2) - r) < eps:
                analytic_boundary_x.append(xi)
                analytic_boundary_y.append(yi)

        return np.array(analytic_boundary_x), np.array(analytic_boundary_y)
