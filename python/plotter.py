from typing import List, Tuple
from pathlib import Path
import math
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation


class Plotter():
    """Parent class of specialized plotters."""
    def __init__(self, folder_path: str,
                 data_drawn: str="phase") -> None:
        self.folder = folder_path
        calc_folder = Path(self.folder) / "calculations"
        if not calc_folder.exists():
            raise FileNotFoundError(f"Složka {calc_folder} neexistuje!")
        self.file_paths = sorted([f for f in calc_folder.iterdir()])

        self.fig = None
        self.ax = None
        self.data_drawn = data_drawn
        self.output_name_tag = None
        self.configuration = None
        config_path = Path(self.folder) / "info" / "config.json"
        self._load_configuration(config_path)

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

    def _load_configuration(self, file_path: Path):
        if not file_path.exists():
            raise FileNotFoundError(f"Soubor {file_path} neexistuje!")

        with file_path.open("r", encoding="utf-8") as file:
            self.configuration = json.load(file)

    def check_for_info_folder(self) -> None:
        path_to_visual = Path(self.folder) / "info"
        if not path_to_visual.exists():
            path_to_visual.mkdir(parents=True, exist_ok=True)
            print("Creating folder \"info\" in ", self.folder)

    def update(self, frame: int):
        pass

    def save_animation(self) -> None:
        print("Preparing animation for", str(self.output_name_tag)[1:])
        self.check_for_info_folder()
        num_frames = len(self.file_paths)
        ani = animation.FuncAnimation(self.fig, self.update, frames=num_frames, blit=True)
        output_file = Path(self.folder) / "info" / ("ACE" + self.output_name_tag + ".gif")
        ani.save(output_file, writer='pillow', fps=max(num_frames/10, 10))
        print("Finished the animation. Saved at:", output_file)

    def save_frame(self, frame: int) -> None:
        if frame >= len(self.file_paths):
            print(f"Frame {frame} is unavailable.")
            return
        self.check_for_info_folder()
        self.update(frame)
        output_file = Path(self.folder) / "info" / ("ACE" + self.output_name_tag + f"_{frame:05d}.jpg")
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
        self.title = self.ax.set_title("Čas: 0")
        self.output_name_tag = "_" + self.data_drawn + "_2D_"

        #self.ax.vlines(np.unique(x), x.min(), x.max(), colors='gray', alpha=0.6)
        #self.ax.hlines(np.unique(y), y.min(), y.max(), colors='gray', alpha=0.6)

        self.draw_function = draw_function
        self.draw_num_boundary = draw_num_boundary
        self.draw_analit_boundary = draw_analit_boundary

        if self.draw_function:
            self.function_sc = self.ax.scatter([], [], c=[], cmap='viridis', s=40)
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
        time_step = (self.configuration["solver"]["final_time"] - self.configuration["solver"]["initial_time"]) / self.configuration["solver"]["frame_num"]
        time_decimal_places = math.ceil(math.log(time_step, 0.1))
        self.title.set_text(f"Čas: {time_step*frame:.{time_decimal_places}f}")
        print("Updating frame: ", frame)
        artist = []

        if self.draw_function:
            self.function_sc.set_array(value)
            artist.append(self.function_sc)

        if self.draw_num_boundary:
            num_boundary_x, num_boundary_y = self.find_boundary_points(x, y, value, threshold=0.5)
            num_boundary_sorted_x, num_boundary_sorted_y = self.sort_boundary_points(num_boundary_x, num_boundary_y)
            self.num_bound_sc.set_data(num_boundary_sorted_x, num_boundary_sorted_y)
            artist.append(self.num_bound_sc)

        if self.draw_analit_boundary:
            analytic_boundary_x, analytic_boundary_y = self.get_analytic_solution(frame*time_step)
            analytic_sorted_x, analytic_sorted_y = self.sort_boundary_points(analytic_boundary_x, analytic_boundary_y)
            self.anal_bound_sc.set_data(analytic_sorted_x, analytic_sorted_y)
            artist.append(self.anal_bound_sc)

        artist.append(self.title)
        return artist

    def find_boundary_points(self, x, y, phase_values, threshold=0.5) -> Tuple[np.array, np.array]:
        boundary_x = []
        boundary_y = []
        # Převod na mřížku
        grid_shape = (self.configuration["solver"]["sizeX"],
                      self.configuration["solver"]["sizeY"])
        grid_x = x.reshape(grid_shape)
        grid_y = y.reshape(grid_shape)
        grid_values = phase_values.reshape(grid_shape)

        # Iterace přes buňky
        has_over = False
        has_bellow = False
        offsets = [(0, 0), (0, 1), (1, 0), (1, 1)]
        for i in range(0, grid_values.shape[0]-1):
            for j in range(0, grid_values.shape[1]-1):
                has_over = False
                has_bellow = False
                for offs in offsets:
                    if grid_values[i + offs[0], j + offs[1]] < threshold:
                        has_bellow = True
                    else:
                        has_over = True

                interpol_points = []
                if has_over and has_bellow:
                    if (grid_values[i, j] - threshold)*(grid_values[i, j+1] - threshold) < 0:
                        interpol_x = (grid_x[i, j]
                                      + (grid_x[i, j+1] - grid_x[i, j])
                                      / (grid_values[i, j+1] - grid_values[i, j])
                                      * (threshold - grid_values[i, j]))
                        interpol_points.append((interpol_x, grid_y[i, j]))
                    if (grid_values[i+1, j] - threshold)*(grid_values[i+1, j+1] - threshold) < 0:
                        interpol_x = (grid_x[i+1, j]
                                      + (grid_x[i+1, j+1] - grid_x[i+1, j])
                                      / (grid_values[i+1, j+1] - grid_values[i+1, j])
                                      * (threshold - grid_values[i+1, j]))
                        interpol_points.append((interpol_x, grid_y[i+1, j]))
                    if (grid_values[i, j] - threshold)*(grid_values[i+1, j] - threshold) < 0:
                        interpol_y = (grid_y[i, j]
                                      + (grid_y[i+1, j] - grid_y[i, j])
                                      / (grid_values[i+1, j] - grid_values[i, j])
                                      * (threshold - grid_values[i, j]))
                        interpol_points.append((grid_x[i, j], interpol_y))
                    if (grid_values[i, j+1] - threshold)*(grid_values[i+1, j+1] - threshold) < 0:
                        interpol_y = (grid_y[i, j+1]
                                      + (grid_y[i+1, j+1] - grid_y[i, j+1])
                                      / (grid_values[i+1, j+1] - grid_values[i, j+1])
                                      * (threshold - grid_values[i, j+1]))
                        interpol_points.append((grid_x[i, j+1], interpol_y))

                    if len(interpol_points) != 2:
                        raise RuntimeError(
                            f"Invalid number of interpolation points: {len(interpol_points)} "
                            f"while computing cell: {i}, {j}"
                        )
                    boundary_x.append((interpol_points[0][0]+interpol_points[1][0])/2)
                    boundary_y.append((interpol_points[0][1]+interpol_points[1][1])/2)
        return boundary_x, boundary_y

    def sort_boundary_points(self, boundary_x: np.array, boundary_y: np.array):
        sorted_x = []
        sorted_y = []

        if len(boundary_x) == 0:
            return  np.array(sorted_x), np.array(sorted_y)

        max_distance = 2*math.sqrt(pow((self.configuration["solver"]["domain"]["x_right"]
                                        - self.configuration["solver"]["domain"]["x_left"])
                                       / (self.configuration["solver"]["sizeX"] - 1), 2)
                                   + pow((self.configuration["solver"]["domain"]["y_right"]
                                          - self.configuration["solver"]["domain"]["y_left"])
                                         / (self.configuration["solver"]["sizeY"] - 1), 2))

        def get_distance(point1: Tuple, point2: Tuple) -> float:
            return math.sqrt(math.pow(point1[0] - point2[0],2) + math.pow(point1[1] - point2[1],2))

        picked_mask = np.full_like(boundary_x, 0)
        for start_index, start_point in enumerate(zip(boundary_x, boundary_y)):
            if picked_mask[start_index] == 0:
                sorted_x.append(start_point[0])
                sorted_y.append(start_point[1])
                picked_mask[start_index] = 1

                previous_point = start_point
                closest_point = (0, 0)
                best_distance = max_distance
                best_index = None
                while closest_point is not None:
                    closest_point = None
                    best_distance = max_distance
                    best_index = None
                    for index, point in enumerate(zip(boundary_x, boundary_y)):
                        distance = get_distance(previous_point, point)
                        if picked_mask[index]==0 and distance <= best_distance:
                            closest_point = point
                            best_distance = distance
                            best_index = index

                    if closest_point is not None:
                        sorted_x.append(closest_point[0])
                        sorted_y.append(closest_point[1])
                        picked_mask[best_index] = 1
                        previous_point = closest_point
        #sorted_x.append(sorted_x[0])
        #sorted_y.append(sorted_y[0])
        return np.array(sorted_x), np.array(sorted_y)

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

    def get_analytic_solution(self, t) -> Tuple[np.array, np.array]:
        r0 = (self.configuration["solver"]["domain"]["x_right"] - self.configuration["solver"]["domain"]["x_left"]) / 4
        coef = 1

        # There is no circle.
        if r0**2 + 2*t < 0:
            return np.array([]), np.array([])
        r = np.sqrt(r0**2 + coef*2*t)

        # Locating boundary points
        analytic_boundary_x, analytic_boundary_y = [], []

        offset_x = (self.configuration["solver"]["domain"]["x_left"] + self.configuration["solver"]["domain"]["x_right"])/2
        offset_y = (self.configuration["solver"]["domain"]["y_left"] + self.configuration["solver"]["domain"]["y_right"])/2
        for phi in np.linspace(0, 2*math.pi, num=1000):
            analytic_boundary_x.append(r * math.cos(phi) + offset_x)
            analytic_boundary_y.append(r * math.sin(phi) + offset_y)

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
        self.ax.set_zlim(min(0, value.min()*1.5), value.max()*1.5)
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("y")
        self.ax.set_zlabel("p")
        self.title = self.ax.set_title("Čas: 0")
        self.output_name_tag = "_" + self.data_drawn + "_3D"

        self.function_surf = self.ax.plot_surface(x_grid, y_grid, val_grid, color='white', edgecolors="black", linewidth=0.25)

    def update(self, frame: int):
        x, y, value = self._load_data(self.file_paths[frame])
        x_grid, y_grid = np.meshgrid(np.unique(x), np.unique(y))
        val_grid = value.reshape(x_grid.shape)

        time_step = (self.configuration["solver"]["final_time"] - self.configuration["solver"]["initial_time"]) / 100
        time_decimal_places = math.ceil(math.log(time_step, 0.1))
        self.title.set_text(f"Čas: {time_step*frame:.{time_decimal_places}f}")
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

        self.ax.set_ylim(min(0, value.min()*1.5), value.max()*1.5)
        self.ax.set_aspect('auto')
        self.title = self.ax.set_title("Čas: 0")
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
        time_step = (self.configuration["solver"]["final_time"] - self.configuration["solver"]["initial_time"]) / 100
        time_decimal_places = math.ceil(math.log(time_step, 0.1))
        self.title.set_text(f"Čas: {time_step*frame:.{time_decimal_places}f}")
        print("Updating frame: ", frame)
        artist = []

        f_var, f_val = self.get_cut(x, y, value)

        self.function.set_data(f_var, f_val)
        artist.append(self.function)

        if self.draw_analytic_solution:
            self.ana_sol.set_data(self.get_analytic_solution(frame*time_step, 4))
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
            mask = (x - cut_value <= distance / 1.5) & (x - cut_value >= 0)
            return y[mask], value[mask]

        if self.axis != "x" or self.axis != "y":
            raise ValueError("Wrong axis in get_cut.")

    def get_analytic_solution(self, t, n_max):
        var_space = None
        if self.axis == "x":
            var_space = np.linspace(self.configuration["solver"]["domain"]["x_left"],
                                    self.configuration["solver"]["domain"]["x_right"],
                                    500)
        elif self.axis == "y":
            var_space = np.linspace(self.configuration["solver"]["domain"]["y_left"],
                                    self.configuration["solver"]["domain"]["y_right"],
                                    500)
        else:
            raise ValueError("Wrong axis given: ", self.axis)

        func_space = np.zeros_like(var_space)
        F = [0, 1/2, 0, 1/4, 0]
        var_left = var_space.min()
        var_right = var_space.max()

        for n in range(1, n_max + 1):
            lambda_n = pow(n * math.pi / (var_right - var_left), 2)
            C_n = pow(1/2, n)
            for index, var in enumerate(var_space):
                func_space[index] += C_n * math.exp(-lambda_n*t) * math.sin(math.sqrt(lambda_n)*(var-var_left))
                func_space[index] += F[n] / lambda_n * (1 - math.exp(-lambda_n * t)) * math.sin(math.sqrt(lambda_n)*(var-var_left))
        return var_space, func_space
