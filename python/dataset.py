"""Loading of the computed results, separated from the plotting."""
from pathlib import Path
from typing import List, Optional, Tuple
import json
import numpy as np


def model_type(configuration: dict) -> str:
    """The model of the run, "allen_cahn", "Zr-Nb_Choud2013" or "Zr-H_Ma2002" (model.type in the config)."""
    return configuration["model"]["type"]


def force_term(configuration: dict) -> Tuple[str, float]:
    """The prescribed force term of the allen_cahn model as (type, size), e.g. ("constant", 0.0);
    the size is 0 for "radial", which has none. ("", 0.0) for the other models."""
    if model_type(configuration) != "allen_cahn":
        return "", 0.0
    section = configuration["model"]["allen_cahn"]
    return section["force_term"], section.get("force_size", 0.0)


def anisotropy(configuration: dict) -> Optional[dict]:
    """The anisotropy {"A", "m", "theta_0"} of the model (the keys A, m, theta_0 of model.<type>),
    None if the model has none."""
    section = configuration["model"][model_type(configuration)]
    if "A" not in section:
        return None
    return {key: section[key] for key in ("A", "m", "theta_0")}


class Dataset():
    """Access to the results of one computation stored in a results folder.

    The folder has to contain the subfolder "calculations" with the frames and the file
    "info/setup/config.json" with the configuration.
    A frame is stored in "Result-NNNNN.npy", an array of the shape (2, sizeY, sizeX), phase first
    and concentration second. The coordinates are not stored, they are computed from the configuration.
    The last loaded frame is kept in memory, so several plotters asking for
    the same frame read the file only once.
    """
    def __init__(self, folder_path) -> None:
        self.folder = Path(folder_path)

        calc_folder = self.folder / "calculations"
        if not calc_folder.exists():
            raise FileNotFoundError(f"Folder {calc_folder} does not exist!")
        # Only the frames written by the solver, nothing else from the folder.
        self.file_paths: List[Path] = self._find_frames(calc_folder)
        if len(self.file_paths) == 0:
            raise FileNotFoundError(f"Folder {calc_folder} does not contain any Result-*.npy files!")

        # Subfolders of "info": the solver writes "setup", the Python scripts write the other two.
        self.setup_folder = self.folder / "info" / "setup"
        self.plots_folder = self.folder / "info" / "plots"
        self.boundary_folder = self.folder / "info" / "boundary"

        config_path = self.setup_folder / "config.json"
        if not config_path.exists():
            raise FileNotFoundError(f"File {config_path} does not exist!")
        with config_path.open("r", encoding="utf-8") as file:
            self.configuration = json.load(file)

        self.time_step = self._calculate_time_step()

        self._grid_xy = None    # coordinates of the grid points, computed when first needed
        self._cached_frame = None
        self._cached_data = None
        self._value_ranges = None   # (min, max) of each field over all frames, computed when first needed

    @staticmethod
    def _find_frames(calc_folder: Path) -> List[Path]:
        """One file per frame, ordered by the frame number."""
        return sorted(calc_folder.glob("Result-*.npy"), key=lambda path: path.name)

    def _calculate_time_step(self) -> float:
        """Time between two saved frames."""
        t0 = self.configuration["solver"]["initial_time"]
        tf = self.configuration["solver"]["final_time"]
        n  = max(len(self.file_paths) - 1, self.configuration["solver"]["frame_num"])
        return (tf - t0) / n

    def time(self, frame: int) -> float:
        """Time of the given frame."""
        return self.configuration["solver"]["initial_time"] + self.time_step * frame

    def get_max_frame_number(self) -> int:
        return len(self.file_paths) - 1

    def expected_frame_count(self) -> int:
        """Number of frames the solver writes according to the configuration (frame_num + 1, including the initial one)."""
        return self.configuration["solver"]["frame_num"] + 1

    def is_complete(self) -> bool:
        """False if some frames are missing (e.g. deleted by DELETE_CALCS, or the run was interrupted)."""
        return len(self.file_paths) >= self.expected_frame_count()

    def load(self, frame: int, data_drawn: str="phase") -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Loads one frame.

        Args:
            frame (int): Index of the frame (0 is the first file).
            data_drawn (str): "phase" or "concentration".

        Returns:
            x: "x" coordinates of the data (flat array, x changes fastest).
            y: "y" coordinates of the data (flat array).
            value: Values of the data selected by data_drawn (flat array).
            The arrays are read-only, because they are shared through the cache.
        """
        if data_drawn not in ("phase", "concentration"):
            raise ValueError("Wrong value for data_drawn given.")

        x, y, phase, concentration = self._get_frame_data(frame)
        return x, y, phase if data_drawn == "phase" else concentration

    def value_range(self, data_drawn: str="phase") -> Tuple[float, float]:
        """
        The smallest and the largest finite value of the field over all frames, so that the plots of
        all frames (and the animations) can share one scale whatever the model and its units.
        NaN and inf are skipped (the last frame of a failed run); (0, 0) if no value is finite.
        All frames are read once, at the first call, for both fields.
        """
        if data_drawn not in ("phase", "concentration"):
            raise ValueError("Wrong value for data_drawn given.")
        if self._value_ranges is None:
            lows = [np.inf, np.inf]
            highs = [-np.inf, -np.inf]
            for path in self.file_paths:
                fields = np.load(path)
                for k in range(2):
                    finite = fields[k][np.isfinite(fields[k])]
                    if finite.size:
                        lows[k] = min(lows[k], float(finite.min()))
                        highs[k] = max(highs[k], float(finite.max()))
            self._value_ranges = {name: (lows[k], highs[k]) if lows[k] <= highs[k] else (0.0, 0.0)
                                  for k, name in enumerate(("phase", "concentration"))}
        return self._value_ranges[data_drawn]

    def _get_frame_data(self, frame: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Returns x, y, phase and concentration of the frame, reading the file only if it is not cached."""
        if self._cached_frame != frame:
            path = self.file_paths[frame]
            x, y = self._grid_coordinates()
            fields = np.load(path)
            expected = (2, self._size_y(), self._size_x())
            if fields.shape != expected:
                raise ValueError(f"File {path} has the shape {fields.shape}, the configuration says {expected}.")
            phase, concentration = fields[0].ravel(), fields[1].ravel()
            for array in (x, y, phase, concentration):
                array.flags.writeable = False
            self._cached_data = (x, y, phase, concentration)
            self._cached_frame = frame
        return self._cached_data

    def _size_x(self) -> int:
        return self.configuration["solver"]["sizeX"]

    def _size_y(self) -> int:
        return self.configuration["solver"]["sizeY"]

    def _grid_coordinates(self) -> Tuple[np.ndarray, np.ndarray]:
        """Flat arrays x and y of all the grid points (x changes fastest), as the solver computes them."""
        if self._grid_xy is None:
            domain = self.configuration["solver"]["domain"]
            size_x, size_y = self._size_x(), self._size_y()
            hx = (domain["x_right"] - domain["x_left"]) / (size_x - 1)
            hy = (domain["y_right"] - domain["y_left"]) / (size_y - 1)
            x_lines = domain["x_left"] + np.arange(size_x) * hx
            y_lines = domain["y_left"] + np.arange(size_y) * hy
            x, y = (a.ravel() for a in np.meshgrid(x_lines, y_lines))
            x.flags.writeable = False
            y.flags.writeable = False
            self._grid_xy = (x, y)
        return self._grid_xy
