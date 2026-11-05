from typing import List, Tuple
from contextlib import ExitStack
from pathlib import Path
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from dataset import Dataset
from boundary import Boundary, join_curves


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

        # Choose the precision according to the step (at most 6 decimal places)
        prec = max(0, -int(math.floor(math.log10(self.time_step)))) if self.time_step != 0 else 0
        prec = min(6, prec)

        # Compute the maximal number of integer digits in the range, for a fixed width
        max_abs = max(abs(t0), abs(tf))
        int_digits = 1 if max_abs < 1 else int(math.floor(math.log10(max_abs))) + 1
        width = int_digits + (1 if min(t0, tf) < 0 else 0) + 1 + prec  # sign + decimal point

        # Monospace makes every character the same width
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

        # The curves of the boundary (of the phase) come from here, the plotter only draws them.
        self.boundary = Boundary(dataset)

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
            cax = divider.append_axes("right", size="5%", pad=0.05)  # same height as the axes
            self.colorbar = self.fig.colorbar(self.function_im, cax=cax, label=r"$c$")
            self.output_name_tag += "f"
        if self.draw_num_boundary:
            self.num_bound_sc, = self.ax.plot([], [], color='orange', linestyle="-", label='Numerical solution')
            self.output_name_tag = self.output_name_tag + "n"
        if self.draw_analit_boundary:
            self.anal_bound_sc, = self.ax.plot([], [], color='blue', linestyle=":", label='Analytic solution')
            self.output_name_tag = self.output_name_tag + "a"

        #self.ax.legend()

    def update(self, frame: int):
        _, _, value = self._load_data(frame)
        self.set_time_title(frame, self.title)

        print("Updating frame: ", frame)
        artist = []

        if self.draw_function:
            self.function_im.set_data(value.reshape(self.grid_shape))
            artist.append(self.function_im)

        if self.draw_num_boundary:
            curves = self.boundary.numerical_curves(frame)
            self.num_bound_sc.set_data(*join_curves(curves, separate_by_nan=True))
            artist.append(self.num_bound_sc)

        if self.draw_analit_boundary:
            curves = self.boundary.analytic_curves(frame)
            self.anal_bound_sc.set_data(*join_curves(curves, separate_by_nan=True))
            artist.append(self.anal_bound_sc)

        artist.append(self.title)
        return artist

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
    """Plotter of a cut through the data: the values along one grid line of the domain."""
    def __init__(self,
                 dataset: Dataset,
                 data_drawn: str="phase",
                 axis: str="x",
                 part: float=0.5):
        """
        Args:
            dataset (Dataset): The computed results.
            data_drawn (str): "phase" or "concentration".
            axis (str): The axis along which the values are plotted. "x": the values along x on a line
                y = const, "y": the values along y on a line x = const.
            part (float): Position of the line as the fraction (from 0 to 1) of the domain in the other
                coordinate. The first grid line at or after this position is used.
        """
        super().__init__(dataset,
                         data_drawn)
        if axis not in ("x", "y"):
            raise ValueError(f"Invalid axis given: {axis}")
        if not 0 <= part <= 1:
            raise ValueError(f"The part has to be between 0 and 1, got {part}.")
        self.axis = axis
        self.part = part

        x, y, value = self._load_data(0)
        # In the data x changes fastest, so the values reshape to (sizeY, sizeX).
        self.grid_shape = (self.configuration["solver"]["sizeY"], self.configuration["solver"]["sizeX"])
        grid_x_lines = np.array(x.reshape(self.grid_shape)[0])
        grid_y_lines = np.array(y.reshape(self.grid_shape)[:, 0])
        if self.axis == "x":
            self.variable_lines, fixed_lines = grid_x_lines, grid_y_lines
            fixed_name = "y"
        else:
            self.variable_lines, fixed_lines = grid_y_lines, grid_x_lines
            fixed_name = "x"

        # The cut is the first grid line that is not before the wanted position (the last one at most).
        wanted_position = fixed_lines[0] + (fixed_lines[-1] - fixed_lines[0]) * self.part
        self.cut_index = min(int(np.searchsorted(fixed_lines, wanted_position, side="left")),
                             len(fixed_lines) - 1)
        self.cut_position = fixed_lines[self.cut_index]

        self.fig, self.ax = plt.subplots()

        self.ax.set_ylim(0, max(0.1, value.max()*1.2))
        self.ax.set_xlim(self.variable_lines[0], self.variable_lines[-1])
        self.ax.set_box_aspect(1/(self.variable_lines[-1] - self.variable_lines[0]))
        self.ax.set_xlabel(rf"${self.axis}$ (${fixed_name} = {self.cut_position:.4g}$)")
        self.ax.set_ylabel(r"$\phi$" if self.data_drawn == "phase" else r"$c$")

        #self.ax.set_aspect('auto')
        self.title = self.ax.set_title("t = 0")
        self.output_name_tag = "_" + self.data_drawn + "_" + self.axis + "-cut"

        self.function, = self.ax.plot(self.variable_lines,
                                      self.get_cut(value),
                                      color='orange',
                                      linestyle='-')

    def update(self, frame: int):
        _, _, value = self._load_data(frame)
        #self.ax.set_ylim(min(0, value.min()*1.5), math.floor(value.max()*1.2 / 0.01)*0.01)
        self.ax.set_ylim(0, max(0.11, value.max()*1.1))

        self.set_time_title(frame, self.title)

        print("Updating frame: ", frame)

        self.function.set_data(self.variable_lines, self.get_cut(value))
        return [self.function, self.title]

    def get_cut(self, value: np.ndarray) -> np.ndarray:
        """The values (a flat array of the whole frame) on the line of the cut, in the order of the increasing coordinate."""
        grid_values = value.reshape(self.grid_shape)
        if self.axis == "x":
            return grid_values[self.cut_index, :]
        return grid_values[:, self.cut_index]


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
    max_fps = 30  # upper limit because of the file size
    num_frames = len(dataset.file_paths)

    # subsample the frames so that the resulting fps does not exceed max_fps
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
