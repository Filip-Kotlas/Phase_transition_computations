"""Loading of the computed results, separated from the plotting."""
from pathlib import Path
from typing import Dict, List, Tuple
import json
import numpy as np


class Dataset():
    """Access to the results of one computation stored in a results folder.

    The folder has to contain the subfolder "calculations" with the frames and the file
    "info/setup/config.json" with the configuration (older results have it in "info/config.json"). A frame is stored in one of two formats:
      - "Result-NNNNN.npy": array of the shape (2, sizeY, sizeX), phase first and concentration second.
        The coordinates are not stored, they are computed from the configuration.
      - "Result-NNNNN.txt": text file with the lines "x y phase concentration" (the older format).
    If both files exist for a frame, the .npy file is used.
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
            raise FileNotFoundError(f"Folder {calc_folder} does not contain any Result-*.npy or Result-*.txt files!")

        # Subfolders of "info": the solver writes "setup", the Python scripts write the other two.
        self.setup_folder = self.folder / "info" / "setup"
        self.plots_folder = self.folder / "info" / "plots"
        self.boundary_folder = self.folder / "info" / "boundary"

        config_path = self.setup_folder / "config.json"
        if not config_path.exists():
            # Results computed before the subfolders existed have the config directly in "info".
            old_config_path = self.folder / "info" / "config.json"
            if not old_config_path.exists():
                raise FileNotFoundError(f"File {config_path} does not exist!")
            config_path = old_config_path
        with config_path.open("r", encoding="utf-8") as file:
            self.configuration = json.load(file)

        self.time_step = self._calculate_time_step()

        self._grid_xy = None    # coordinates for the .npy frames, computed when first needed
        self._cached_frame = None
        self._cached_data = None

    @staticmethod
    def _find_frames(calc_folder: Path) -> List[Path]:
        """One file per frame, ordered by the frame number. The .npy file wins over the .txt one."""
        by_name: Dict[str, Path] = {path.stem: path for path in calc_folder.glob("Result-*.txt")}
        by_name.update({path.stem: path for path in calc_folder.glob("Result-*.npy")})
        return [by_name[name] for name in sorted(by_name)]

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
        if data_drawn == "phase":
            column = 2
        elif data_drawn == "concentration":
            column = 3
        else:
            raise ValueError("Wrong value for data_drawn given.")

        x, y, phase, concentration = self._get_frame_data(frame)
        return x, y, phase if data_drawn == "phase" else concentration

    def _get_frame_data(self, frame: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Returns x, y, phase and concentration of the frame, reading the file only if it is not cached."""
        if self._cached_frame != frame:
            path = self.file_paths[frame]
            if path.suffix == ".npy":
                x, y = self._grid_coordinates()
                fields = np.load(path)
                expected = (2, self._size_y(), self._size_x())
                if fields.shape != expected:
                    raise ValueError(f"File {path} has the shape {fields.shape}, the configuration says {expected}.")
                phase, concentration = fields[0].ravel(), fields[1].ravel()
            else:
                data = np.loadtxt(path)
                x, y, phase, concentration = data[:, 0], data[:, 1], data[:, 2], data[:, 3]
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
