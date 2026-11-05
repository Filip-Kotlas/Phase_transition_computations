from typing import List, Tuple
from contextlib import ExitStack
from pathlib import Path
import math
import numpy as np
import contourpy
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

        # Coordinates of the grid lines (in the data x changes fastest). They are copied, so the whole
        # frame is not kept in memory by them.
        grid_shape = (self.configuration["solver"]["sizeY"], self.configuration["solver"]["sizeX"])
        self.grid_x_lines = np.array(x.reshape(grid_shape)[0])
        self.grid_y_lines = np.array(y.reshape(grid_shape)[:, 0])

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
            num_curves = self.find_boundary_curves(x, y, value, threshold=0.5)
            self.num_bound_sc.set_data(*self.join_curves(num_curves, separate_by_nan=True))
            artist.append(self.num_bound_sc)
            self.save_boundary_points("interface_num.crv", *self.join_curves(num_curves), frame)
            #Carefull the boundary is saved during animation and also during showing saving one frame. Could pose problems.
            # TODO: solve that

        if self.draw_analit_boundary:
            analytic_curves = self.get_analytic_curves(frame*self.time_step)
            self.anal_bound_sc.set_data(*self.join_curves(analytic_curves, separate_by_nan=True))
            artist.append(self.anal_bound_sc)
            self.save_boundary_points("interface_anal.crv", *self.join_curves(analytic_curves), frame)

        artist.append(self.title)
        return artist

    def find_boundary_curves(self, x, y, phase_values, threshold=0.5) -> List[np.ndarray]:
        """
        Finds the boundary, i.e. the level "threshold" of the phase, as a list of curves.

        The phase is interpolated linearly along the edges of the grid. The points of the curves are the points
        on the edges where the interpolated phase equals the threshold. The curves are found by contourpy,
        which returns the points already ordered along every curve. A closed curve ends with its first point
        again. A curve that reaches the border of the domain is open.

        Returns:
            List of arrays of the shape (number of points, 2), the columns are x and y.
        """
        # In the data x changes fastest, so the shape is (sizeY, sizeX): the first index
        # of the grid is the position in y, the second one the position in x.
        grid_shape = (self.configuration["solver"]["sizeY"],
                      self.configuration["solver"]["sizeX"])
        grid_x = x.reshape(grid_shape)
        grid_y = y.reshape(grid_shape)
        grid_values = phase_values.reshape(grid_shape)

        generator = contourpy.contour_generator(grid_x[0], grid_y[:, 0], grid_values)
        return generator.lines(threshold)

    def join_curves(self, curves: List[np.ndarray], separate_by_nan: bool=False) -> Tuple[np.ndarray, np.ndarray]:
        """
        Puts the curves one after another into one pair of arrays (x, y).

        A closed curve already ends with its first point. With separate_by_nan a NaN point is inserted
        between the curves, so a line that is plotted from the result is interrupted there.
        """
        parts = []
        for curve in curves:
            if separate_by_nan and len(parts) > 0:
                parts.append(np.full((1, 2), np.nan))
            parts.append(curve)
        if len(parts) == 0:
            return np.array([]), np.array([])
        joined = np.vstack(parts)
        return joined[:, 0], joined[:, 1]
    
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

    def get_analytic_curves(self, t) -> List[np.ndarray]:
        """
        The boundary of the analytic solution at the time t as a list of curves, like find_boundary_curves.

        The points of the curves are the exact intersections of the analytic boundary with the lines of the grid
        (x = const and y = const through the grid points). So they lie on the same edges of the grid as the points
        of the numerical boundary and they are spaced in the same way. Only the part of the boundary inside
        the domain is returned. The list is empty if there is no shape.
        """
        smallest_step = min(np.diff(self.grid_x_lines).min(), np.diff(self.grid_y_lines).min())

        # The exact boundary is replaced by a polyline whose segments are much shorter than the grid step,
        # so the error of this replacement is negligible.
        number_of_segments = 2000
        while True:
            curve_x, curve_y = self.get_analytic_solution(t, number_of_segments)
            if len(curve_x) == 0:
                return []
            if np.hypot(np.diff(curve_x), np.diff(curve_y)).max() <= smallest_step / 50:
                break
            number_of_segments *= 2
            if number_of_segments > 100_000_000:
                raise RuntimeError("The analytic boundary cannot be approximated by a fine enough polyline.")
        return self.crossings_with_grid_lines(curve_x, curve_y)

    def crossings_with_grid_lines(self, curve_x, curve_y) -> List[np.ndarray]:
        """
        Finds the points where a closed polyline crosses the lines of the grid, ordered along the polyline.

        The polyline has to be closed (its last point is its first point) and its segments have to be shorter
        than the grid step. Only the points inside the domain are kept. If the whole polyline is inside
        the domain, the result is one closed curve (it ends with its first point). Otherwise the result are open
        curves, one for every part of the polyline that lies inside the domain.
        """
        grid_x, grid_y = self.grid_x_lines, self.grid_y_lines
        n = len(curve_x) - 1  # number of the segments

        # Number of the grid lines that are not above the coordinate, i.e. between which lines the coordinate is.
        strip_x = np.searchsorted(grid_x, curve_x, side="right")
        strip_y = np.searchsorted(grid_y, curve_y, side="right")
        if np.abs(np.diff(strip_x)).max() > 1 or np.abs(np.diff(strip_y)).max() > 1:
            raise RuntimeError("A segment of the polyline crosses more lines of the grid, it is too long.")

        # Crossings with the vertical lines (x = const): the segments whose ends lie in different strips.
        segment_v = np.nonzero(strip_x[:-1] != strip_x[1:])[0]
        line_x = grid_x[np.minimum(strip_x[segment_v], strip_x[segment_v + 1])]
        position_v = (line_x - curve_x[segment_v]) / (curve_x[segment_v + 1] - curve_x[segment_v])
        point_v = np.column_stack((line_x, curve_y[segment_v]
                                   + position_v * (curve_y[segment_v + 1] - curve_y[segment_v])))
        # Crossings with the horizontal lines (y = const).
        segment_h = np.nonzero(strip_y[:-1] != strip_y[1:])[0]
        line_y = grid_y[np.minimum(strip_y[segment_h], strip_y[segment_h + 1])]
        position_h = (line_y - curve_y[segment_h]) / (curve_y[segment_h + 1] - curve_y[segment_h])
        point_h = np.column_stack((curve_x[segment_h]
                                   + position_h * (curve_x[segment_h + 1] - curve_x[segment_h]), line_y))

        # All the crossings in the order along the polyline: by the segment, inside a segment by the position.
        segment = np.concatenate((segment_v, segment_h))
        points = np.vstack((point_v, point_h))
        order = np.lexsort((np.concatenate((position_v, position_h)), segment))
        segment, points = segment[order], points[order]

        def in_domain(x, y):
            return (x >= grid_x[0]) & (x <= grid_x[-1]) & (y >= grid_y[0]) & (y <= grid_y[-1])

        inside = in_domain(curve_x, curve_y)  # for the points of the polyline
        if not inside.any():
            return []

        closed = inside.all()
        if closed:
            groups = [points]
        else:
            # The parts of the polyline inside the domain are numbered. A crossing belongs to the part of
            # its segment that is inside, crossings outside of the domain are dropped.
            inside_vertex = inside[:n]
            run = np.cumsum(inside_vertex & ~np.roll(inside_vertex, 1))
            first_inside, second_inside = inside[segment], inside[segment + 1]
            point_run = np.where(first_inside, run[segment], run[(segment + 1) % n])
            keep = (first_inside | second_inside) & in_domain(points[:, 0], points[:, 1])
            points, point_run = points[keep], point_run[keep]

            groups = [points[point_run == number] for number in range(1, run.max() + 1)]
            if inside_vertex[0] and inside_vertex[-1]:
                # The last part continues in the first one (the numbering starts at its beginning from 0).
                groups[-1] = np.vstack((groups[-1], points[point_run == 0]))

        curves = []
        for group in groups:
            # The same point can be found twice (the polyline goes exactly through a point of the grid).
            different = np.ones(len(group), dtype=bool)
            different[1:] = np.any(group[1:] != group[:-1], axis=1)
            group = group[different]
            if closed:
                if len(group) > 1 and np.all(group[-1] == group[0]):
                    group = group[:-1]
                if len(group) > 0:
                    group = np.vstack((group, group[:1]))  # closed by the first point
            if len(group) >= 2:
                curves.append(group)
        return curves

    def get_analytic_solution(self, t, number_of_segments=1000) -> Tuple[np.array, np.array]:
        """The analytic boundary as a closed polyline with number_of_segments segments (x and y coordinates)."""
        # The radius is relative to the shorter side of the domain, the same as in the solver (Parameters.cu).
        domain = self.configuration["solver"]["domain"]
        smaller_side = min(domain["x_right"] - domain["x_left"],
                           domain["y_right"] - domain["y_left"])
        r0 = self.configuration["init_cond"]["radius_proportion"] * smaller_side
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

        offset_x = (self.configuration["solver"]["domain"]["x_left"] + self.configuration["solver"]["domain"]["x_right"])/2
        offset_y = (self.configuration["solver"]["domain"]["y_left"] + self.configuration["solver"]["domain"]["y_right"])/2
        phi = np.linspace(0, 2*math.pi, num=number_of_segments, endpoint=False)
        analytic_boundary_x = r * np.cos(phi) + offset_x
        analytic_boundary_y = r * np.sin(phi) + offset_y

        # The curve is closed by its first point.
        return (np.append(analytic_boundary_x, analytic_boundary_x[0]),
                np.append(analytic_boundary_y, analytic_boundary_y[0]))

    def get_analytic_solution_wulff_shape(self, t, r0, number_of_segments=1000) -> Tuple[np.array, np.array]:
        def psi(theta) -> float:
            theta_0 = self.configuration["anisotropy"]["theta_0"]
            m = self.configuration["anisotropy"]["m"]
            A = self.configuration["anisotropy"]["A"]
            return 1 + A*np.cos(m*(theta - theta_0))
        
        def der_psi(theta) -> float:
            theta_0 = self.configuration["anisotropy"]["theta_0"]
            m = self.configuration["anisotropy"]["m"]
            A = self.configuration["anisotropy"]["A"]
            return -m*A*np.sin(m*(theta - theta_0))
        
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

        offset_x = (self.configuration["solver"]["domain"]["x_left"] + self.configuration["solver"]["domain"]["x_right"])/2
        offset_y = (self.configuration["solver"]["domain"]["y_left"] + self.configuration["solver"]["domain"]["y_right"])/2
        theta = np.linspace(0, 2*math.pi, num=number_of_segments, endpoint=False)
        analytic_boundary_x = r*(psi(theta)*np.cos(theta) - der_psi(theta)*np.sin(theta)) + offset_x
        analytic_boundary_y = r*(psi(theta)*np.sin(theta) + der_psi(theta)*np.cos(theta)) + offset_y

        # The curve is closed by its first point.
        return (np.append(analytic_boundary_x, analytic_boundary_x[0]),
                np.append(analytic_boundary_y, analytic_boundary_y[0]))

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
