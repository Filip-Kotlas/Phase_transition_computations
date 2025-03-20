from typing import List, Tuple
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


class Plotter():
    """Parent class of specialized plotters."""
    def __init__(self, folder_path: str,
                 data_drawn: str="phase") -> None:
        self.folder = folder_path
        self.file_paths = sorted([os.path.join(self.folder, "calculations", f)
                                  for f
                                  in os.listdir(os.path.join(folder_path, "calculations"))
                                  if f.endswith('.txt')])
        self.fig = None
        self.ax = None
        self.data_drawn = data_drawn
        self.output_name_tag = None

    def _load_data(self, file_path: str) -> Tuple[List[float], List[float], List[float]]:
        """
        Loads data from a file.
        
        Args:
            file_path (str): Path of the file with data.

        Returns:
            x: List of "x" coordinates of the data.
            y: List of "y" coordinates of the data.
            value: List of values of the data based on the self.data_drawn variable.
        """
        data = np.loadtxt(file_path)
        x = data[:, 0]
        y = data[:, 1]
        if self.data_drawn == "phase":
            value = data[:, 2]
        elif self.data_drawn == "concentration":
            value = data[:, 3]
        else:
            raise ValueError("Wrong value for data_drawn given.")
        return x, y, value

    def check_for_info_folder(self) -> None:
        path_to_visual = os.path.join(self.folder, "info")
        if not os.path.exists(path_to_visual):
            os.makedirs(path_to_visual)
            print("Creating folder \"info\" in ", self.folder)

    def update(self, frame: int):
        pass

    def save_animation(self) -> None:
        self.check_for_info_folder()
        num_frames = len(self.file_paths)
        ani = animation.FuncAnimation(self.fig, self.update, frames=num_frames, blit=True)
        output_file = os.path.join(self.folder, "info", "ACE" + self.output_name_tag + ".gif")
        ani.save(output_file, writer='pillow', fps=max(num_frames/10, 10))
        print("Finished the animation. Saved at:", output_file)

    def save_frame(self, frame: int) -> None:
        if frame >= len(self.file_paths):
            print(f"Frame {frame} is unavailable.")
            return
        self.check_for_info_folder()
        self.update(frame)
        output_file = os.path.join(self.folder, "info", "ACE" + self.output_name_tag + f"_{frame:05d}.jpg")
        self.fig.savefig(output_file, dpi=300)
        print(f"Frame {frame} saved at:", output_file)

    def show_frame(self, frame: int) -> None:
        if frame >= len(self.file_paths):
            print(f"Frame {frame} is unavailable.")
            return
        self.update(frame)
        #self.ax.axis("off")
        plt.show()

class BoundaryPlotter2D(Plotter):
    """Plotter of a 2D plot."""
    def __init__(self,
                 folder_path: str,
                 draw_function: bool=False,
                 draw_num_boundary: bool=True,
                 draw_analit_boundary: bool=False,
                 data_drawn: str="phase"):
        super().__init__(folder_path,
                         data_drawn)

        x, y, _ = self._load_data(self.file_paths[0])
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(x.min(), x.max())
        self.ax.set_ylim(y.min(), y.max())
        self.ax.set_aspect('equal')
        self.title = self.ax.set_title("Frame: 0")
        self.output_name_tag = "_" + self.data_drawn + "_2D_"

        self.draw_function = draw_function
        self.draw_num_boundary = draw_num_boundary
        self.draw_analit_boundary = draw_analit_boundary

        if self.draw_function:
            self.function_sc = self.ax.scatter([], [], c=[], cmap='viridis', s=1)
            self.function_sc.set_offsets(np.column_stack((x, y)))
            self.output_name_tag = self.output_name_tag + "f"
        if self.draw_num_boundary:
            self.num_bound_sc, = self.ax.plot([], [], color='orange', linestyle="-", label='Num. sol.')
            self.output_name_tag = self.output_name_tag + "n"
        if self.draw_analit_boundary:
            self.anal_bound_sc, = self.ax.plot([], [], color='blue', linestyle=":", label='Anal. sol.')
            self.output_name_tag = self.output_name_tag + "a"

        self.ax.legend()

    def update(self, frame: int):
        x, y, value = self._load_data(self.file_paths[frame])
        self.title.set_text(f"Frame: {frame}")
        print("Updating frame: ", frame)
        artist = []

        if self.draw_function:
            self.function_sc.set_array(value)
            artist.append(self.function_sc)

        if self.draw_num_boundary:
            num_boundary_x, num_boundary_y = self.find_boundary(x, y, value, threshold=0.5)
            self.num_bound_sc.set_data(num_boundary_x, num_boundary_y)
            artist.append(self.num_bound_sc)

        if self.draw_analit_boundary:
            analytic_boundary_x, analytic_boundary_y = self.get_analytic_solution(x, y, 0.5, frame*0.001, (x[1]-x[0])/2)
            self.anal_bound_sc.set_data(analytic_boundary_x, analytic_boundary_y)
            artist.append(self.anal_bound_sc)

        artist.append(self.title)
        return artist

    def find_boundary(self, x, y, phase_values, threshold=0.5) -> Tuple[np.array, np.array]:
        # Předpokládáme, že body jsou na mřížce (nebo téměř na mřížce)
        grid_shape = (len(np.unique(x)), len(np.unique(y)))

        # Převod na matici (2D mřížku)
        grid_x = x.reshape(grid_shape)
        grid_y = y.reshape(grid_shape)
        grid_values = phase_values.reshape(grid_shape)

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

    def reset_draw_settings(self,
                            draw_function: bool=False,
                            draw_num_boundary: bool=True,
                            draw_anal_boundary: bool=False,
                            data_drawn: str="phase") -> None:
        self.draw_function = draw_function
        self.draw_num_boundary = draw_num_boundary
        self.draw_analit_boundary = draw_anal_boundary

        self.data_drawn = data_drawn
        self.output_name_tag = "_" + self.data_drawn + "_3D_"
        if self.draw_function:
            self.output_name_tag = self.output_name_tag + "f"
        if self.draw_num_boundary:
            self.output_name_tag = self.output_name_tag + "n"
        if self.draw_analit_boundary:
            self.output_name_tag = self.output_name_tag + "a"

    def get_analytic_solution(self, x, y, r0, t, eps) -> Tuple[np.array, np.array]:
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

class SurfacePlotter(Plotter):
    """Plotter of graph of a function of two variables."""

    def __init__(self,
                 folder_path,
                 data_drawn: str="concentration"):
        super().__init__(folder_path, data_drawn)

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')

        x, y, value = self._load_data(self.file_paths[0])
        x_grid, y_grid = np.meshgrid(np.unique(x), np.unique(y))
        val_grid = value.reshape(x_grid.shape)

        self.ax.set_xlim(x.min(), x.max())
        self.ax.set_ylim(y.min(), y.max())
        self.ax.set_zlim(value.min()*1.5, value.max()*1.5)
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_zlabel("p")
        self.title = self.ax.set_title("Frame: 0")
        self.output_name_tag = "_" + self.data_drawn + "_3D"

        self.function_surf = self.ax.plot_surface(x_grid, y_grid, val_grid, color='white', edgecolors="black", linewidth=0.25)

    def update(self, frame: int):
        x, y, value = self._load_data(self.file_paths[frame])
        x_grid, y_grid = np.meshgrid(np.unique(x), np.unique(y))
        val_grid = value.reshape(x_grid.shape)

        self.title.set_text(f"Frame: {frame}")
        print("Updating frame: ", frame)
        artist = []

        alpha = 1

        self.function_surf.remove()
        self.function_surf = self.ax.plot_surface(x_grid,
                                                  y_grid,
                                                  val_grid,
                                                  color='white',
                                                  edgecolors="black",
                                                  linewidth=0.25,
                                                  alpha=alpha)

        artist.append(self.title)
        return artist

class CutPlotter(Plotter):
    """Plotter of a 2D scatter plot."""
    def __init__(self,
                 folder_path: str,
                 data_drawn: str="phase",
                 axis: str="x",
                 part: float=0.5,
                 draw_analytic_solution: bool=False):
        super().__init__(folder_path,
                         data_drawn)
        self.axis = axis
        self.part = part
        self.draw_analytic_solution = draw_analytic_solution

        x, y, value = self._load_data(self.file_paths[0])
        self.fig, self.ax = plt.subplots()

        if self.axis == "x":
            self.ax.set_xlim(x.min(), x.max())
        elif self.axis == "y":
            self.ax.set_xlim(y.min(), y.max())
        else:
            raise ValueError("Invalid axis given")

        self.ax.set_ylim(value.min()*1.5, value.max()*1.5)
        self.ax.set_aspect('equal')
        self.title = self.ax.set_title("Frame: 0")
        self.output_name_tag = "_" + self.data_drawn + "_" + self.axis + "-cut"

        self.function, = self.ax.plot([],
                                      [],
                                      color='orange',
                                      linestyle='-',
                                      label=self.axis + "-cut of c at " + str(int(self.part*100)) + "%")
        if self.axis == "x":
            self.function.set_data(x, value)
        elif self.axis == "y":
            self.function.set_data(y, value)

        if self.draw_analytic_solution:
            self.ana_sol, = self.ax.plot([],
                                         [],
                                        color='b',
                                        linestyle=':',
                                        label="Analytic solution")

        self.ax.legend()

    def update(self, frame: int):
        x, y, value = self._load_data(self.file_paths[frame])
        self.title.set_text(f"Frame: {frame}")
        print("Updating frame: ", frame)
        artist = []

        f_var, f_val = self.get_cut(x, y, value)

        self.function.set_data(f_var, f_val)
        artist.append(self.function)

        if self.draw_analytic_solution:
            self.ana_sol.set_data(f_var, self.get_analytic_solution(frame*0.001, f_var, 5))
            artist.append(self.ana_sol)

        artist.append(self.title)
        return artist

    def get_cut(self, x, y, value):
        if self.axis == "x":
            cut_value = y.min() + (y.max() - y.min()) * self.part
            unique_y = np.unique(y)
            distance = unique_y[1] - unique_y[0]
            mask = (y - cut_value <= distance / 1.5) & (y - cut_value >= 0)
            return x[mask], value[mask]

        if self.axis == "y":
            cut_value = x.min() + (x.max() - x.min()) * self.part
            unique_x = np.unique(x)
            distance = unique_x[1] - unique_x[0]
            mask = abs(x - cut_value) < distance / 2
            return y[mask], value[mask]

        if self.axis != "x" or self.axis != "y":
            raise ValueError("Wrong axis in get_cut.")

    def get_analytic_solution(self, t, x_arr, n):
        result = np.zeros(len(x_arr))
        for k in range(1, n):
            lambda_k = pow(k * math.pi * 2 / (x_arr.max() - x_arr.min()), 2)
            C_k = pow(1/2, k)
            for index, x in enumerate(x_arr):
                result[index] += C_k * math.exp(-lambda_k*t) * math.sin(math.sqrt(lambda_k)*x)
        return result
