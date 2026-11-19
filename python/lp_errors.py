"""L2 and L-infinity errors of a run against a reference run (usually a finer one), frame by frame.

Usage:
    python3 python/lp_errors.py RUN REFERENCE [--eval-cells NX NY] [--output FILE]

The Python version of lp_dist_new_new.pas from Measurement/ (M. Benes). For every frame and for
both fields (phase and concentration):
    1. Both runs are interpolated bilinearly (the functions Inter1 and Inter2 of the Pascal code),
       so that runs on different grids can be compared: I_run and I_ref.
    2. The difference d = I_run - I_ref is evaluated at the points of the evaluation grid ("fine
       mesh" in the Pascal code), which covers the same domain:
           L-infinity error = max |d| over the points of the evaluation grid,
           L2 error         = sqrt( sum w_ij * d_ij^2 * hx_eval * hy_eval ), the composite
                              trapezoidal rule (weight 1 inside, 1/2 on the edges, 1/4 in the corners).
Then the errors of all frames are combined into norms in time:
    L-infinity in time = max over the frames (what the Pascal code writes into max.txt),
    L2 in time         = sqrt( integral over (initial_time, final_time) of error(t)^2 dt ), the
                         trapezoidal rule over the frames (not in the Pascal code).

The evaluation grid is by default the least common multiple of the cell counts of both runs, in
each direction (e.g. 800 cells for runs with 400 and 800 cells). It contains the nodes of both
grids, which makes the L-infinity error exact (d is bilinear on each of its cells, so |d| has its
maximum at a node). --eval-cells sets another grid, e.g. the 2000 x 2000 used with the Pascal code.

The runs have to cover the same domain, have the same initial and final time, the same number of
frames and the same model, and both have to be finished with all their frames. Any other difference
in the settings the solver uses (physical parameters, boundary or initial conditions, ...) only gives
a warning; parts of the configs the solver does not use are not compared (see used_settings).
Only the current config layout (with the section "model") is supported; results with the old layout
can be converted by python3 python/migrate_config.py (it changes their config file in place).
The concentration is not measured for models that do not compute it (allen_cahn).

Lines of the report that are not data rows start with "#", so the file written by --output can be
read e.g. by numpy.loadtxt or gnuplot. Exit code: 0 on success, 1 when the runs cannot be compared,
2 for wrong arguments.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import argparse
import math
import sys
import numpy as np
from dataset import Dataset
from compare_runs import flatten, read_setup

FIELDS = ("phase", "concentration")
# The config sections the solver reads (Parameters::load), besides model.type and model.<type>.
SOLVER_SECTIONS = ("solver", "init_cond", "bound_cond")
# Keys that are compared by check_compatible() itself (they have to be equal) or that differ on purpose.
NOT_WARNED_KEYS = ("solver.sizeX", "solver.sizeY")
# The same list as in cli.py (importing cli.py would load matplotlib and the plotters).
MODELS_WITHOUT_CONCENTRATION = ("allen_cahn",)
# The default evaluation grid (the least common multiple of the cell counts) may be at most this
# many times finer than the finer of the two runs; otherwise --eval-cells has to be given.
MAX_EVALUATION_REFINEMENT = 4
# Number of evaluation points interpolated at once. It limits the memory: 16 MB per float64 array of
# a block, about ten such arrays at once (the corner values, t_1, t_2, the difference, the weights).
POINTS_PER_BLOCK = 2_000_000


class RunError(Exception):
    """The runs cannot be compared (missing, unfinished, incompatible, NaN in a frame)."""


@dataclass
class FieldErrors:
    """Errors of one field, one value per frame."""
    l2: np.ndarray
    linf: np.ndarray


@dataclass
class ErrorResult:
    run: Path
    reference: Path
    model: str
    run_cells: Tuple[int, int]          # (cells in x, cells in y) = (sizeX - 1, sizeY - 1)
    reference_cells: Tuple[int, int]
    eval_cells: Tuple[int, int]
    eval_grid_is_default: bool
    h: float                            # max(hx, hy) of the run, the h of the EOC (as in the Pascal code)
    times: np.ndarray                   # time of every frame
    time_step: float                    # time between two frames
    fields: Dict[str, FieldErrors]      # only the measured fields
    warnings: List[str] = field(default_factory=list)


def max_in_time(errors: np.ndarray) -> float:
    """L-infinity norm in time: the largest error over the frames."""
    return float(np.max(errors))


def l2_in_time(errors: np.ndarray, time_step: float) -> float:
    """L2 norm in time: sqrt of the integral of error(t)^2, trapezoidal rule over the frames."""
    weights = np.ones(len(errors))
    weights[0] = weights[-1] = 0.5
    return float(np.sqrt(time_step * np.sum(weights * errors**2)))


class BilinearInterpolation:
    """Bilinear interpolation of the values of one grid at the points of the evaluation grid.

    The same formula as Inter1/Inter2 in lp_dist_new_new.pas. The coordinates are measured from
    the lower left corner of the domain. For a point (x, y):
        i0 = trunc(x/h_x), s = x/h_x - i0         (the cell containing the point and the position in it)
        j0 = trunc(y/h_y), r = y/h_y - j0
        T_1 = (T[i0+1, j0  ] - T[i0, j0  ])*s + T[i0, j0  ]    (along x at the bottom of the cell)
        T_2 = (T[i0+1, j0+1] - T[i0, j0+1])*s + T[i0, j0+1]    (along x at the top of the cell)
        value = T_1 + (T_2 - T_1)*r
    T[i, j] is the value at the node (i, j); in our arrays of the shape (sizeY, sizeX) it is values[j, i].
    On the right or top edge of the domain i0 (or j0) is the last node and s (or r) is 0. The Pascal
    code then reads the node behind the edge, which is 0 there, multiplied by 0; the values are padded
    with a row and a column of zeros here, so the result is the same, bit for bit.
    """
    def __init__(self, cells: Tuple[int, int], h_x: float, h_y: float,
                 eval_x: np.ndarray, eval_y: np.ndarray) -> None:
        self.cells = cells
        x_scaled = eval_x / h_x
        y_scaled = eval_y / h_y
        self.i0 = np.trunc(x_scaled).astype(np.int64)
        self.j0 = np.trunc(y_scaled).astype(np.int64)
        if self.i0.min() < 0 or self.i0.max() > cells[0] or self.j0.min() < 0 or self.j0.max() > cells[1]:
            raise RunError("an evaluation point lies outside the grid (internal error)")
        self.s = x_scaled - self.i0
        self.r = y_scaled - self.j0

    def pad(self, values: np.ndarray) -> np.ndarray:
        """The values of the shape (sizeY, sizeX) with a row and a column of zeros added behind the edges."""
        padded = np.zeros((self.cells[1] + 2, self.cells[0] + 2))
        padded[:-1, :-1] = values
        return padded

    def evaluate(self, padded: np.ndarray, rows: slice) -> np.ndarray:
        """The interpolated values (padded by pad()) at the evaluation points of the given rows (y),
        an array of the shape (number of rows, all points in x)."""
        i0 = self.i0[np.newaxis, :]
        j0 = self.j0[rows, np.newaxis]
        s = self.s[np.newaxis, :]
        r = self.r[rows, np.newaxis]
        t_1 = (padded[j0, i0 + 1] - padded[j0, i0]) * s + padded[j0, i0]
        t_2 = (padded[j0 + 1, i0 + 1] - padded[j0 + 1, i0]) * s + padded[j0 + 1, i0]
        return t_1 + (t_2 - t_1) * r


def cells_of(dataset: Dataset) -> Tuple[int, int]:
    solver = dataset.configuration["solver"]
    return solver["sizeX"] - 1, solver["sizeY"] - 1


def steps_of(dataset: Dataset) -> Tuple[float, float]:
    """The space steps (hx, hy) of the run (h1, h2 of the Pascal code)."""
    domain = dataset.configuration["solver"]["domain"]
    cells = cells_of(dataset)
    return (domain["x_right"] - domain["x_left"]) / cells[0], (domain["y_right"] - domain["y_left"]) / cells[1]


def open_run(path: Path, label: str) -> Dataset:
    """Loads the run and checks that it finished with all its frames."""
    if not (path / "calculations").is_dir():
        raise RunError(f"the {label} {path} does not exist or has no calculations/ folder")
    try:
        dataset = Dataset(path)
    except (FileNotFoundError, ValueError) as error:
        raise RunError(f"the {label} {path}: {error}")
    if "model" not in dataset.configuration:
        raise RunError(f"the {label} {path} has the old config layout (no section 'model'), which is not supported; "
                       "python3 python/migrate_config.py can convert its config file (in place)")
    # A run that recorded its start but not its end crashed, was killed, failed or is still running.
    setup = read_setup(path)
    if setup.has_run_date and setup.computation_seconds is None:
        raise RunError(f"the {label} {path} did not finish (no 'Computation time' in parameters.txt)")
    if not dataset.is_complete():
        raise RunError(f"the {label} {path} has {dataset.get_max_frame_number() + 1} frames, "
                       f"the config says {dataset.expected_frame_count()}")
    return dataset


def used_settings(configuration: dict) -> Dict[str, object]:
    """The settings the solver uses, as flattened config keys with their values ({"solver.sizeX": 200, ...}).

    The solver reads the sections solver, init_cond and bound_cond, model.type and the section of the
    chosen model, model.<type>; everything else in the config (the sections of the other models, the
    old-layout sections kept in tests/configs, ...) is not used. Of the read keys, two groups are
    required by the solver but not used:
        solver.integration_time_step, when solver.custom_integration_time_step is false
            (the solver starts from min(hx, hy)^2 / 5),
        init_cond.init_cond_file_path and init_cond.scaling, when init_cond.init_cond_from_file is false.
    """
    model = configuration["model"]["type"]
    used: Dict[str, object] = {"model.type": model}
    for section in SOLVER_SECTIONS:
        used.update(flatten(configuration.get(section, {}), section + "."))
    used.update(flatten(configuration["model"].get(model, {}), f"model.{model}."))
    if used.get("solver.custom_integration_time_step") is False:
        used.pop("solver.integration_time_step", None)
    if used.get("init_cond.init_cond_from_file") is False:
        used.pop("init_cond.init_cond_file_path", None)
        used.pop("init_cond.scaling", None)
    return used


def check_compatible(run: Dataset, reference: Dataset) -> List[str]:
    """Raises RunError if the runs cannot be compared; returns the warnings about other differences."""
    run_solver, ref_solver = run.configuration["solver"], reference.configuration["solver"]
    problems = []
    for key in ("x_left", "x_right", "y_left", "y_right"):
        if run_solver["domain"][key] != ref_solver["domain"][key]:
            problems.append(f"domain.{key} {run_solver['domain'][key]} vs {ref_solver['domain'][key]}")
    for key in ("initial_time", "final_time", "frame_num"):
        if run_solver[key] != ref_solver[key]:
            problems.append(f"{key} {run_solver[key]} vs {ref_solver[key]}")
    run_model, ref_model = run.configuration["model"]["type"], reference.configuration["model"]["type"]
    if run_model != ref_model:
        problems.append(f"model {run_model} vs {ref_model}")
    if problems:
        raise RunError("the runs differ (run vs reference): " + ", ".join(problems))

    # A setting used in one run only (e.g. integration_time_step, when custom_integration_time_step
    # is true in one run only) is a difference too.
    warnings = []
    run_used, ref_used = used_settings(run.configuration), used_settings(reference.configuration)
    for key in sorted(set(run_used) | set(ref_used)):
        if key in NOT_WARNED_KEYS:
            continue
        run_value = run_used.get(key, "not used")
        ref_value = ref_used.get(key, "not used")
        if run_value != ref_value:
            warnings.append(f"config differs: {key} = {run_value} (run) vs {ref_value} (reference)")
    run_cells, ref_cells = cells_of(run), cells_of(reference)
    if ref_cells[0] <= run_cells[0] or ref_cells[1] <= run_cells[1]:
        warnings.append(f"the reference ({ref_cells[0]} x {ref_cells[1]} cells) is not finer than "
                        f"the run ({run_cells[0]} x {run_cells[1]} cells)")
    return warnings


def evaluation_cells(run_cells: Tuple[int, int], ref_cells: Tuple[int, int],
                     requested: Optional[Sequence[int]]) -> Tuple[Tuple[int, int], List[str]]:
    """The cells of the evaluation grid and the warnings about it."""
    if requested is not None:
        cells = (int(requested[0]), int(requested[1]))
        if min(cells) < 1:
            raise RunError("--eval-cells has to be positive")
        warnings = []
        for axis, name in ((0, "x"), (1, "y")):
            if cells[axis] % run_cells[axis] != 0 or cells[axis] % ref_cells[axis] != 0:
                warnings.append(f"the evaluation grid ({cells[axis]} cells in {name}) does not contain all nodes of "
                                f"both grids ({run_cells[axis]} and {ref_cells[axis]} cells), so the "
                                "L-infinity error can be underestimated")
        return cells, warnings
    cells = (math.lcm(run_cells[0], ref_cells[0]), math.lcm(run_cells[1], ref_cells[1]))
    for axis, name in ((0, "x"), (1, "y")):
        finest = max(run_cells[axis], ref_cells[axis])
        if cells[axis] > MAX_EVALUATION_REFINEMENT * finest:
            raise RunError(f"the least common multiple of the cell counts in {name} ({run_cells[axis]} and "
                           f"{ref_cells[axis]}) is {cells[axis]}, more than {MAX_EVALUATION_REFINEMENT}x the finer grid; "
                           "choose the evaluation grid with --eval-cells NX NY")
    return cells, []


def frame_fields(dataset: Dataset, frame: int, label: str) -> Dict[str, np.ndarray]:
    """Phase and concentration of a frame as float64 arrays of the shape (sizeY, sizeX)."""
    size_x, size_y = dataset.configuration["solver"]["sizeX"], dataset.configuration["solver"]["sizeY"]
    result = {}
    for name in FIELDS:
        try:
            _, _, values = dataset.load(frame, name)
        except ValueError as error:
            raise RunError(f"the {label} {dataset.folder}: {error}")
        values = values.astype(np.float64).reshape(size_y, size_x)
        if not np.all(np.isfinite(values)):
            raise RunError(f"frame {frame} of the {label} {dataset.folder} contains NaN or inf")
        result[name] = values
    return result


def compute_errors(run_path: Path, reference_path: Path,
                   eval_cells: Optional[Sequence[int]] = None) -> ErrorResult:
    """The errors of the run against the reference in every frame. Raises RunError."""
    run = open_run(run_path, "run")
    reference = open_run(reference_path, "reference")
    warnings = check_compatible(run, reference)

    run_cells, ref_cells = cells_of(run), cells_of(reference)
    cells, grid_warnings = evaluation_cells(run_cells, ref_cells, eval_cells)
    warnings += grid_warnings

    # The steps as in the Pascal code: h1 is the step of the run, the steps of the reference
    # (h1_ = N1*h1/N1_) and of the evaluation grid (h1__ = N1*h1/N1__) follow from the domain N1*h1.
    run_hx, run_hy = steps_of(run)
    ref_hx, ref_hy = run_cells[0] * run_hx / ref_cells[0], run_cells[1] * run_hy / ref_cells[1]
    eval_hx, eval_hy = run_cells[0] * run_hx / cells[0], run_cells[1] * run_hy / cells[1]
    eval_x = np.arange(cells[0] + 1) * eval_hx
    eval_y = np.arange(cells[1] + 1) * eval_hy
    run_interpolation = BilinearInterpolation(run_cells, run_hx, run_hy, eval_x, eval_y)
    ref_interpolation = BilinearInterpolation(ref_cells, ref_hx, ref_hy, eval_x, eval_y)

    # Trapezoidal weights: 1 inside, 1/2 on the edges, 1/4 in the corners (= weight_x * weight_y).
    weight_x = np.ones(cells[0] + 1)
    weight_x[0] = weight_x[-1] = 0.5
    weight_y = np.ones(cells[1] + 1)
    weight_y[0] = weight_y[-1] = 0.5

    model = run.configuration["model"]["type"]
    measured = [name for name in FIELDS
                if not (name == "concentration" and model in MODELS_WITHOUT_CONCENTRATION)]
    solver = run.configuration["solver"]
    frame_count = solver["frame_num"] + 1
    time_step = (solver["final_time"] - solver["initial_time"]) / solver["frame_num"]
    times = solver["initial_time"] + time_step * np.arange(frame_count)

    errors = {name: FieldErrors(np.zeros(frame_count), np.zeros(frame_count)) for name in measured}
    rows_per_block = max(1, POINTS_PER_BLOCK // (cells[0] + 1))
    for frame in range(frame_count):
        run_values = frame_fields(run, frame, "run")
        ref_values = frame_fields(reference, frame, "reference")
        for name in measured:
            run_padded = run_interpolation.pad(run_values[name])
            ref_padded = ref_interpolation.pad(ref_values[name])
            linf, sum_squares = 0.0, 0.0
            for start in range(0, cells[1] + 1, rows_per_block):
                rows = slice(start, min(start + rows_per_block, cells[1] + 1))
                difference = (run_interpolation.evaluate(run_padded, rows)
                              - ref_interpolation.evaluate(ref_padded, rows))
                linf = max(linf, float(np.max(np.abs(difference))))
                weights = weight_y[rows, np.newaxis] * weight_x[np.newaxis, :]
                sum_squares += float(np.sum(weights * difference**2))
            errors[name].linf[frame] = linf
            errors[name].l2[frame] = math.sqrt(sum_squares * eval_hx * eval_hy)

    return ErrorResult(run=run_path, reference=reference_path, model=model,
                       run_cells=run_cells, reference_cells=ref_cells,
                       eval_cells=cells, eval_grid_is_default=eval_cells is None,
                       h=max(run_hx, run_hy), times=times, time_step=time_step,
                       fields=errors, warnings=warnings)


def format_report(result: ErrorResult) -> str:
    """The report: a header, one row per frame and the norms in time. All lines but the rows start with #."""
    def number(value: float) -> str:
        return f"{value:.8e}"

    lines = [f"# Run:       {result.run}   ({result.run_cells[0]} x {result.run_cells[1]} cells, h = {result.h:.6e})",
             f"# Reference: {result.reference}   ({result.reference_cells[0]} x {result.reference_cells[1]} cells)",
             f"# Evaluation grid: {result.eval_cells[0]} x {result.eval_cells[1]} cells"
             + (" (the least common multiple of the cell counts)" if result.eval_grid_is_default else " (--eval-cells)"),
             f"# Model: {result.model}"]
    if "concentration" not in result.fields:
        lines.append(f"# The model {result.model} does not compute the concentration, only the phase is measured.")
    lines += [f"# Warning: {warning}" for warning in result.warnings]

    labels = {"phase": "phase", "concentration": "conc"}
    columns = ["frame", "time"]
    for name in result.fields:
        columns += [f"L2_{labels[name]}", f"Linf_{labels[name]}"]
    widths = [6, 16] + [16] * (len(columns) - 2)
    lines.append("#" + "".join(f"{column:>{width}}" for column, width in zip(columns, widths))[1:])
    for frame, time in enumerate(result.times):
        values = [f"{frame:>6}", f"{number(time):>16}"]
        for errors in result.fields.values():
            values += [f"{number(errors.l2[frame]):>16}", f"{number(errors.linf[frame]):>16}"]
        lines.append("".join(values))

    lines.append("# Norms in time (L2 in time: trapezoidal rule over the frames):")
    for name, errors in result.fields.items():
        lines.append(f"#   {name}:")
        lines.append(f"#     L-infinity(0,T; L2)         = {number(max_in_time(errors.l2))}")
        lines.append(f"#     L-infinity(0,T; L-infinity) = {number(max_in_time(errors.linf))}")
        lines.append(f"#     L2(0,T; L2)                 = {number(l2_in_time(errors.l2, result.time_step))}")
        lines.append(f"#     L2(0,T; L-infinity)         = {number(l2_in_time(errors.linf, result.time_step))}")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="L2 and L-infinity errors of a run against a reference run, frame by frame "
                    "(the Python version of lp_dist_new_new.pas).")
    parser.add_argument("run", type=Path, help="the run whose error is measured")
    parser.add_argument("reference", type=Path, help="the reference run, usually a finer one")
    parser.add_argument("--eval-cells", type=int, nargs=2, metavar=("NX", "NY"),
                        help="cells of the evaluation grid (default: the least common multiple of the cell counts)")
    parser.add_argument("--output", type=Path, help="also write the report into this file")
    args = parser.parse_args()

    try:
        result = compute_errors(args.run, args.reference, args.eval_cells)
    except RunError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    report = format_report(result)
    print(report)
    if args.output is not None:
        args.output.write_text(report + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
