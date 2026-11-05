"""Loading of the computed results, separated from the plotting."""
from pathlib import Path
from typing import List, Tuple
import json
import numpy as np


class Dataset():
    """Access to the results of one computation stored in a results folder.

    The folder has to contain the subfolder "calculations" with the frames
    "Result-NNNNN.txt" and the file "info/config.json" with the configuration.
    The last loaded frame is kept in memory, so several plotters asking for
    the same frame read the file only once.
    """
    def __init__(self, folder_path) -> None:
        self.folder = Path(folder_path)

        calc_folder = self.folder / "calculations"
        if not calc_folder.exists():
            raise FileNotFoundError(f"Složka {calc_folder} neexistuje!")
        # Only the frames written by the solver, nothing else from the folder.
        self.file_paths: List[Path] = sorted(calc_folder.glob("Result-*.txt"))
        if len(self.file_paths) == 0:
            raise FileNotFoundError(f"Složka {calc_folder} neobsahuje žádné soubory Result-*.txt!")

        config_path = self.folder / "info" / "config.json"
        if not config_path.exists():
            raise FileNotFoundError(f"Soubor {config_path} neexistuje!")
        with config_path.open("r", encoding="utf-8") as file:
            self.configuration = json.load(file)

        self.time_step = self._calculate_time_step()

        self._cached_frame = None
        self._cached_data = None

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

        data = self._get_frame_data(frame)
        return data[:, 0], data[:, 1], data[:, column]

    def _get_frame_data(self, frame: int) -> np.ndarray:
        """Returns all four columns of the frame, reading the file only if it is not cached."""
        if self._cached_frame != frame:
            data = np.loadtxt(self.file_paths[frame])
            data.flags.writeable = False
            self._cached_data = data
            self._cached_frame = frame
        return self._cached_data
