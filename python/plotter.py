from typing import List, Tuple
from contextlib import ExitStack
from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from dataset import Dataset


class Plotter():
    """Parent class of specialized plotters."""
    def __init__(self, dataset: Dataset,
                 data_drawn: str="phase") -> None:
        self.dataset = dataset
        self.folder = self.dataset.folder
        self.file_paths = self.dataset.file_paths
        self.configuration = self.dataset.configuration

        self.fig = None
        self.ax = None
        self.data_drawn = data_drawn
        self.output_name_tag = None

        self.time_step = self.dataset.time_step

    def _load_data(self, frame: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Loads data of one frame.

        Args:
            frame (int): Index of the frame.

        Returns:
            x: "x" coordinates of the data.
            y: "y" coordinates of the data.
            value: Values of the data based on the self.data_drawn variable.
            The arrays are read-only (they are shared through the cache of the dataset).
        """
        return self.dataset.load(frame, self.data_drawn)

    def get_max_frame_number(self) -> int:
        return self.dataset.get_max_frame_number()

    def check_for_info_folder(self) -> None:
        path_to_visual = Path(self.folder) / "info"
        if not path_to_visual.exists():
            path_to_visual.mkdir(parents=True, exist_ok=True)
            print("Creating folder \"info\" in ", self.folder)

    def set_time_title(self, frame: int, title) -> None:
        t0 = self.configuration["solver"]["initial_time"]
        tf = self.configuration["solver"]["final_time"]
        t = self.dataset.time(frame)

        # Zvol přesnost podle kroku (max 6 desetin)
        prec = max(0, -int(math.floor(math.log10(self.time_step)))) if self.time_step != 0 else 0
        prec = min(6, prec)

        # Spočítej maximální počet celých číslic v rozsahu pro fixní šířku
        max_abs = max(abs(t0), abs(tf))
        int_digits = 1 if max_abs < 1 else int(math.floor(math.log10(max_abs))) + 1
        width = int_digits + (1 if min(t0, tf) < 0 else 0) + 1 + prec  # znaménko + tečka

        # Monospace zajistí, že každý znak má stejnou šířku
        title.set_fontfamily('monospace')
        title.set_text(f"t = {t:>{width}.{prec}f}")

    def update(self, frame: int):
        pass

    def save_animation(self, file_format: str="auto") -> None:
        save_animations([self], file_format)

    def save_frame(self, frame: int) -> None:
        if frame >= len(self.file_paths):
            print(f"Frame {frame} is unavailable.")
            return
        self.check_for_info_folder()
        self.update(frame)
        output_file = Path(self.folder) / "info" / ("ACE" + self.output_name_tag + f"_{frame:05d}.jpg")
        self.fig.savefig(output_file, dpi=300, bbox_inches='tight')
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
                 dataset: Dataset,
                 draw_function: bool=False,
                 draw_num_boundary: bool=True,
                 draw_analit_boundary: bool=False,
                 data_drawn: str="phase"):
        super().__init__(dataset,
                         data_drawn)

        x, y, _ = self._load_data(0)
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlabel(r"$x$")
        self.ax.set_ylabel(r"$y$")
        self.ax.set_xlim(x.min(), x.max())
        self.ax.set_ylim(y.min(), y.max())
        self.ax.set_aspect('equal')
        self.title = self.ax.set_title("t = 0")
        self.output_name_tag = "_" + self.data_drawn + "_2D_"

        #self.ax.vlines(np.unique(x), x.min(), x.max(), colors='gray', alpha=0.6)
        #self.ax.hlines(np.unique(y), y.min(), y.max(), colors='gray', alpha=0.6)

        self.draw_function = draw_function
        self.draw_num_boundary = draw_num_boundary
        self.draw_analit_boundary = draw_analit_boundary

        if self.draw_function:
            # The values are drawn as an image, one pixel per grid point.
            # In the data x changes fastest, so the values reshape to (sizeY, sizeX).
            n_x = self.configuration["solver"]["sizeX"]
            n_y = self.configuration["solver"]["sizeY"]
            self.grid_shape = (n_y, n_x)
            # The extent is given by the edges of the image, the grid points are the
            # centres of the pixels, hence the half of a cell is added on each side.
            half_hx = (x.max() - x.min()) / (n_x - 1) / 2
            half_hy = (y.max() - y.min()) / (n_y - 1) / 2
            extent = (x.min() - half_hx, x.max() + half_hx,
                      y.min() - half_hy, y.max() + half_hy)
            self.function_im = self.ax.imshow(np.zeros(self.grid_shape), origin="lower", extent=extent,
                                              cmap='viridis', vmin=0.0, vmax=0.09)
            divider = make_axes_locatable(self.ax)
            cax = divider.append_axes("right", size="5%", pad=0.05)  # stejné výšky
            self.colorbar = self.fig.colorbar(self.function_im, cax=cax, label=r"$c$")
            self.output_name_tag += "f"
        if self.draw_num_boundary:
            self.num_bound_sc, = self.ax.plot([], [], color='orange', linestyle="-", label='Numerické řešení')
            self.output_name_tag = self.output_name_tag + "n"
        if self.draw_analit_boundary:
            self.anal_bound_sc, = self.ax.plot([], [], color='blue', linestyle=":", label='Analytické řešení')
            self.output_name_tag = self.output_name_tag + "a"

        #self.ax.legend()

    def update(self, frame: int):
        x, y, value = self._load_data(frame)
        self.set_time_title(frame, self.title)

        print("Updating frame: ", frame)
        artist = []

        if self.draw_function:
            self.function_im.set_data(value.reshape(self.grid_shape))
            artist.append(self.function_im)

        if self.draw_num_boundary:
            num_boundary_x, num_boundary_y = self.find_boundary_points(x, y, value, threshold=0.5)
            num_boundary_sorted_x, num_boundary_sorted_y = self.sort_boundary_points(num_boundary_x, num_boundary_y)
            self.num_bound_sc.set_data(num_boundary_sorted_x, num_boundary_sorted_y)
            artist.append(self.num_bound_sc)
            self.save_boundary_points("interface_num.crv", num_boundary_sorted_x, num_boundary_sorted_y, frame)
            #Carefull the boundary is saved during animation and also during showing saving one frame. Could pose problems.
            # TODO: solve that

        if self.draw_analit_boundary:
            analytic_boundary_x, analytic_boundary_y = self.get_analytic_solution(frame*self.time_step)
            analytic_sorted_x, analytic_sorted_y = self.sort_boundary_points(analytic_boundary_x, analytic_boundary_y)
            self.anal_bound_sc.set_data(analytic_sorted_x, analytic_sorted_y)
            artist.append(self.anal_bound_sc)
            self.save_boundary_points("interface_anal.crv", analytic_sorted_x, analytic_sorted_y, frame)

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
        sorted_x.append(sorted_x[0])
        sorted_y.append(sorted_y[0])
        return np.array(sorted_x), np.array(sorted_y)
    
    def save_boundary_points(
        self,
        file_name: str,
        boundary_x: np.array,
        boundary_y: np.array,
        frame: int,
    ) -> None:
        self.check_for_info_folder()
        output_file = Path(self.folder) / "info" / file_name
        time_level = self.dataset.time(frame)

        mode = "w" if frame == 0 else "a"
        with output_file.open(mode, encoding="utf-8") as file:
            file.write(f"# Time level T =  {time_level:.6f}\n")
            for x, y in zip(boundary_x, boundary_y):
                file.write(f"    {x:.17g}   {y:.17g}\n")
            file.write("##\n\n")

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
        r0 = (self.configuration["init_cond"]["radius_proportion"] 
             * (self.configuration["solver"]["domain"]["x_right"] - self.configuration["solver"]["domain"]["x_left"]))
        number_of_segments = 1000
        if self.configuration["init_cond"]["initial_condition"] == "ConstantCircle":
            return self.get_analytic_solution_circle(t, r0, number_of_segments)
        elif self.configuration["init_cond"]["initial_condition"] == "WulffShape":
            return self.get_analytic_solution_wulff_shape(t, r0, number_of_segments)
        else:
            raise ValueError("Wrong initial condition for analytic solution.")

    def get_analytic_solution_circle(self, t, r0, number_of_segments=1000) -> Tuple[np.array, np.array]:
        def r_squared_for_force_0(t) -> float:
            return r0**2 - 2*t

        def r_squared_for_force_radial(t) -> float:
            return r0**2 + 2*t

        def r_used(t) -> float:
            pass

        if self.configuration["force_term"]["type"] == "constant" and self.configuration["force_term"]["size"] == 0:
            r_used = r_squared_for_force_0
        elif self.configuration["force_term"]["type"] == "radial":
            r_used = r_squared_for_force_radial
        else:
            return np.array([]), np.array([])

        # There is no shape.
        if r_used(t) < 0:
            return np.array([]), np.array([])
        r = math.sqrt(r_used(t))

        # Locating boundary points
        analytic_boundary_x, analytic_boundary_y = [], []

        offset_x = (self.configuration["solver"]["domain"]["x_left"] + self.configuration["solver"]["domain"]["x_right"])/2
        offset_y = (self.configuration["solver"]["domain"]["y_left"] + self.configuration["solver"]["domain"]["y_right"])/2
        for phi in np.linspace(0, 2*math.pi, num=number_of_segments):
            analytic_boundary_x.append(r * math.cos(phi) + offset_x)
            analytic_boundary_y.append(r * math.sin(phi) + offset_y)

        analytic_boundary_x.append(analytic_boundary_x[0])
        analytic_boundary_y.append(analytic_boundary_y[0])

        return np.array(analytic_boundary_x), np.array(analytic_boundary_y)

    def get_analytic_solution_wulff_shape(self, t, r0, number_of_segments=1000) -> Tuple[np.array, np.array]:
        def psi(theta) -> float:
            theta_0 = self.configuration["anisotropy"]["theta_0"]
            m = self.configuration["anisotropy"]["m"]
            A = self.configuration["anisotropy"]["A"]
            return 1 + A*math.cos(m*(theta - theta_0))
        
        def der_psi(theta) -> float:
            theta_0 = self.configuration["anisotropy"]["theta_0"]
            m = self.configuration["anisotropy"]["m"]
            A = self.configuration["anisotropy"]["A"]
            return -m*A*math.sin(m*(theta - theta_0))
        
        def r_squared_for_force_0(t) -> float:
            return r0**2 - 2*t

        def r_squared_for_force_radial(t) -> float:
            return r0**2 + 2*t

        def r_used(t) -> float:
            pass

        if self.configuration["force_term"]["type"] == "constant" and self.configuration["force_term"]["size"] == 0:
            r_used = r_squared_for_force_0
        elif self.configuration["force_term"]["type"] == "radial":
            r_used = r_squared_for_force_radial
        else:
            return np.array([]), np.array([])

        # There is no shape.
        if r_used(t) < 0:
            return np.array([]), np.array([])
        r = math.sqrt(r_used(t))

        theta_0 = self.configuration["anisotropy"]["theta_0"]
        # Locating boundary points
        analytic_boundary_x, analytic_boundary_y = [], []

        offset_x = (self.configuration["solver"]["domain"]["x_left"] + self.configuration["solver"]["domain"]["x_right"])/2
        offset_y = (self.configuration["solver"]["domain"]["y_left"] + self.configuration["solver"]["domain"]["y_right"])/2
        for theta in np.linspace(0, 2*math.pi, num=number_of_segments):
            x = r*(psi(theta)*math.cos(theta) - der_psi(theta)*math.sin(theta))
            y = r*(psi(theta)*math.sin(theta) + der_psi(theta)*math.cos(theta))
            analytic_boundary_x.append(x + offset_x)
            analytic_boundary_y.append(y + offset_y)

        analytic_boundary_x.append(analytic_boundary_x[0])
        analytic_boundary_y.append(analytic_boundary_y[0])

        return np.array(analytic_boundary_x), np.array(analytic_boundary_y)

class SurfacePlotter(Plotter):
    """Plotter of graph of a function of two variables."""
    def __init__(self,
                 dataset: Dataset,
                 data_drawn: str="concentration"):
        super().__init__(dataset, data_drawn)

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')

        x, y, value = self._load_data(0)
        x_grid, y_grid = np.meshgrid(np.unique(x), np.unique(y))
        val_grid = value.reshape(x_grid.shape)

        self.ax.set_xlim(x.min(), x.max())
        self.ax.set_ylim(y.min(), y.max())
        if( data_drawn == "phase"):
            self.ax.set_zlim(0, 1.1)
            self.ax.set_box_aspect([x.max() - x.min(), y.max() - y.min(), 1])
        elif( data_drawn == "concentration"):
            self.ax.set_zlim(0, 0.09)
            self.ax.set_box_aspect([x.max() - x.min(), y.max() - y.min(), 0.5])
        else:
            raise ValueError("Wrong data_drawn given.")

        self.ax.set_xlabel(r"$x$")
        self.ax.set_ylabel(r"$y$")
        if data_drawn == "concentration":
            self.ax.set_zlabel(r"$c$")
        elif data_drawn == "phase":
            self.ax.set_zlabel(r"$\phi$")
        else:
            raise ValueError("Wrong type of plot data.")

        self.title = self.ax.set_title("t = 0")
        self.output_name_tag = "_" + self.data_drawn + "_3D"

        self.function_surf = self.ax.plot_surface(x_grid, y_grid, val_grid, color='white', edgecolors="black", linewidth=0.25)

    def update(self, frame: int):
        x, y, value = self._load_data(frame)
        x_grid, y_grid = np.meshgrid(np.unique(x), np.unique(y))
        val_grid = value.reshape(x_grid.shape)

        self.set_time_title(frame, self.title)

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
                 dataset: Dataset,
                 data_drawn: str="phase",
                 axis: str="x",
                 part: float=0.5,
                 draw_analytic_solution: bool=False):
        super().__init__(dataset,
                         data_drawn)
        self.axis = axis
        self.part = part
        self.draw_analytic_solution = draw_analytic_solution

        x, y, value = self._load_data(0)
        self.fig, self.ax = plt.subplots()

        self.ax.set_ylim(0, max(0.1, value.max()*1.2))
        if self.axis == "x":
            self.ax.set_xlim(x.min(), x.max())
            self.ax.set_box_aspect(1/(x.max() - x.min()))
        elif self.axis == "y":
            self.ax.set_xlim(y.min(), y.max())
            self.ax.set_box_aspect(1/(y.max() - y.min()))
        else:
            raise ValueError("Invalid axis given")

        #self.ax.set_aspect('auto')
        self.title = self.ax.set_title("t = 0")
        self.output_name_tag = "_" + self.data_drawn + "_" + self.axis + "-cut"

        self.function, = self.ax.plot([],
                                      [],
                                      color='orange',
                                      linestyle='-')
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

        #self.ax.legend()

    def update(self, frame: int):
        x, y, value = self._load_data(frame)
        #self.ax.set_ylim(min(0, value.min()*1.5), math.floor(value.max()*1.2 / 0.01)*0.01)
        self.ax.set_ylim(0, max(0.11, value.max()*1.1))

        self.set_time_title(frame, self.title)

        print("Updating frame: ", frame)
        artist = []

        f_var, f_val = self.get_cut(x, y, value)

        self.function.set_data(f_var, f_val)
        artist.append(self.function)

        if self.draw_analytic_solution:
            self.ana_sol.set_data(self.get_analytic_solution(frame*self.time_step, 4))
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


def ffmpeg_available() -> bool:
    """
    Tells whether matplotlib is able to write mp4 files, i.e. whether ffmpeg can be found.

    The ffmpeg from the PATH is preferred. If there is none, the ffmpeg bundled in the python
    package imageio-ffmpeg is used (when the package is installed).
    """
    if animation.FFMpegWriter.isAvailable():
        return True
    try:
        import imageio_ffmpeg
    except ImportError:
        return False
    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
    return animation.FFMpegWriter.isAvailable()


def save_animations(plotters: List[Plotter], file_format: str="auto") -> None:
    """
    Saves the animations of all the given plotters at once.

    The animations are created frame by frame and in every frame all the plotters are updated.
    Therefore the file with the data of the frame is read only once for all the plotters,
    not once for every plotter (it is so, because the plotters share one Dataset).

    Args:
        plotters (List[Plotter]): Plotters whose animations are saved. They have to share the same Dataset.
        file_format (str): "mp4", "gif" or "auto". With "auto" the mp4 is used when ffmpeg is available,
            otherwise the gif. The mp4 is much smaller and its frames are not collected in memory.
    """
    if file_format not in ("auto", "mp4", "gif"):
        raise ValueError("Wrong file_format given, use auto, mp4 or gif.")
    if len(plotters) == 0:
        return
    dataset = plotters[0].dataset
    if any(plotter.dataset is not dataset for plotter in plotters):
        raise ValueError("All the plotters have to share the same dataset.")

    if file_format == "auto":
        if ffmpeg_available():
            file_format = "mp4"
        else:
            file_format = "gif"
            print("ffmpeg was not found, the animations are saved as gif.")
    elif file_format == "mp4" and not ffmpeg_available():
        raise RuntimeError("ffmpeg was not found, mp4 cannot be saved. Install it or use the gif format.")

    if file_format == "mp4":
        # The codec (h264) needs the width and the height of the frame (in pixels) to be even numbers.
        for plotter in plotters:
            width, height = plotter.fig.canvas.get_width_height()
            if width % 2 != 0 or height % 2 != 0:
                raise ValueError(f"The size of the figure is {width} x {height} pixels, "
                                 "for mp4 both numbers have to be even.")

    duration_s = 10
    max_fps = 30  # horní limit kvůli velikosti souboru
    num_frames = len(dataset.file_paths)

    # podvzorkování snímků tak, aby výsledné fps nepřekročilo max_fps
    step = max(1, int(np.ceil(num_frames / (duration_s * max_fps))))
    frames = list(range(0, num_frames, step))

    fps = max(1, int(round(len(frames) / duration_s)))

    output_files = []
    for plotter in plotters:
        print("Preparing animation for", str(plotter.output_name_tag)[1:])
        plotter.check_for_info_folder()
        output_files.append(Path(plotter.folder) / "info" / ("ACE" + plotter.output_name_tag + "." + file_format))

    # The gif writer collects the frames in memory and creates the file when the "with" block ends.
    # The mp4 writer sends every frame to ffmpeg at once.
    with ExitStack() as stack:
        writers = []
        for plotter, output_file in zip(plotters, output_files):
            if file_format == "mp4":
                # yuv420p: the pixel format that is played by all the players (also PowerPoint, browsers),
                # crf 18: the quality of the compression (lower is better, 18 is visually almost lossless).
                writer = animation.FFMpegWriter(fps=fps, codec="libx264",
                                                extra_args=["-pix_fmt", "yuv420p", "-crf", "18"])
            else:
                writer = animation.PillowWriter(fps=fps)
            stack.enter_context(writer.saving(plotter.fig, output_file, plotter.fig.dpi))
            writers.append(writer)

        for frame in frames:
            for plotter, writer in zip(plotters, writers):
                plotter.update(frame)
                writer.grab_frame()

    for output_file in output_files:
        print("Finished the animation. Saved at:", output_file)
