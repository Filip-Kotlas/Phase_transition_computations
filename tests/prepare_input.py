"""Writes the initial-condition file that a test config reads, if it reads one.

Usage: python3 tests/prepare_input.py CONFIG RUN_FOLDER

If the config has init_cond_from_file: true, the file is written to
RUN_FOLDER/calculations/<init_cond_file_path> (the solver resolves the path relative to that
folder). Otherwise nothing happens. The file is generated rather than stored in the repo, and
it is always the same, so both versions of the code in a regression test read identical input.

The file is a Result-*.npy frame stored as float32. The solver converts float and double frames to
its Real type, so the file works for float and double builds. (Versions of the solver before
2026-11-11 accept only the Real type of the build, i.e. float32 for float builds; a regression
test against such a version as the baseline needs a float build.)
"""
from pathlib import Path
import json
import sys
import numpy as np

# Equilibrium concentrations of the two phases, as constants::c_init_alpha / c_init_beta in include/constants.hpp.
C_ALPHA = 0.007
C_BETA = 0.025


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit("Usage: prepare_input.py CONFIG RUN_FOLDER")
    config_path, run_folder = Path(sys.argv[1]), Path(sys.argv[2])
    config = json.loads(config_path.read_text(encoding="utf-8"))

    init_cond = config["init_cond"]
    if not init_cond.get("init_cond_from_file", False):
        return
    target = run_folder / "calculations" / init_cond["init_cond_file_path"]
    if target.suffix != ".npy":
        sys.exit(f"{config_path}: init_cond_file_path has to be a .npy file, the solver reads no other format.")

    solver = config["solver"]
    domain = solver["domain"]
    size_x, size_y = solver["sizeX"], solver["sizeY"]
    x_lines = np.linspace(domain["x_left"], domain["x_right"], size_x)
    y_lines = np.linspace(domain["y_left"], domain["y_right"], size_y)
    x, y = np.meshgrid(x_lines, y_lines)   # shape (size_y, size_x), x changes fastest

    # A smooth circle (alpha inside) placed off the centre, so that a mistake in the placement
    # of the file data (e.g. a transposed or shifted grid) changes the result.
    width = domain["x_right"] - domain["x_left"]
    height = domain["y_right"] - domain["y_left"]
    centre_x = domain["x_left"] + 0.45 * width
    centre_y = domain["y_left"] + 0.55 * height
    radius = 0.3 * min(width, height)
    distance = np.hypot(x - centre_x, y - centre_y)
    phase = 0.5 * (1.0 - np.tanh((distance - radius) / (0.03 * min(width, height))))
    concentration = C_ALPHA * phase + C_BETA * (1.0 - phase)

    # The Result-*.npy format: shape (2, sizeY, sizeX), phase first; x changes fastest, as in phase[j, i].
    target.parent.mkdir(parents=True, exist_ok=True)
    np.save(target, np.stack([phase, concentration]).astype("<f4"))


if __name__ == "__main__":
    main()
