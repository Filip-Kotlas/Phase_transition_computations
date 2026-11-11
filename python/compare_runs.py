"""Compares the results of two runs frame by frame.

Usage:
    python3 python/compare_runs.py RUN_A RUN_B [--tol 1e-4]

RUN_A and RUN_B are either two run folders (each with calculations/ and info/setup/), e.g.
results/Test_2 and results/Test_3, or two folders whose subfolders are runs. In the second case
the subfolders with the same name are compared, e.g. results/MyTest/base and results/MyTest/new
of a regression test (tests/regression.sh), where each subfolder is one test case.

For each pair of runs it checks that both runs finished and have the same grid and the same
number of frames, and then compares the frames that both have:
    - whether the frame is bitwise identical,
    - max |phase_A - phase_B|,
    - max |conc_A - conc_B| divided by the range of the concentration in RUN_A (over all frames),
    - the number of cells with phase > 0.5 (alpha) in both runs; a change means that the
      interface moved by at least one cell somewhere.
It also prints the computation time and the git commit of both runs (from parameters.txt).

Verdict of a pair:
    IDENTICAL     all frames bitwise identical and the same number of frames
    CLOSE         differences at most --tol (the phase absolutely, the concentration relative to
                  its range) and the same number of alpha cells in every frame
    DIFFERENT     anything else
    FAILED        a run is missing, did not finish, or contains NaN
    INCOMPATIBLE  the grids differ, so the frames cannot be compared

The last line of the output is the overall count of the verdicts (cluster.sh list shows it).
Exit code: 0 when every pair is IDENTICAL or CLOSE, 1 otherwise, 2 for wrong arguments.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import argparse
import re
import sys
import numpy as np
from dataset import Dataset

VERDICTS = ["IDENTICAL", "CLOSE", "DIFFERENT", "FAILED", "INCOMPATIBLE"]
DEFAULT_TOLERANCE = 1e-4


@dataclass
class SetupInfo:
    """What the binary (and make run or tests/regression.sh) wrote into parameters.txt."""
    found: bool = False
    has_run_date: bool = False
    computation_seconds: Optional[float] = None
    commit: Optional[str] = None


@dataclass
class FrameRow:
    frame: int
    identical: bool
    max_dphase: float
    max_dconc: float          # absolute; divided by the concentration range when printed
    alpha_cells_a: int
    alpha_cells_b: int


@dataclass
class PairResult:
    name: str
    path_a: Path
    path_b: Path
    verdict: str = "FAILED"
    notes: List[str] = field(default_factory=list)
    rows: List[FrameRow] = field(default_factory=list)
    conc_scale: float = 1.0
    setup_a: SetupInfo = field(default_factory=SetupInfo)
    setup_b: SetupInfo = field(default_factory=SetupInfo)

    # np.max, not the built-in max: it keeps a NaN, which max() would silently skip.
    def max_dphase(self) -> Optional[float]:
        return float(np.max([r.max_dphase for r in self.rows])) if self.rows else None

    def max_dconc_relative(self) -> Optional[float]:
        return float(np.max([r.max_dconc for r in self.rows])) / self.conc_scale if self.rows else None


def is_run(path: Path) -> bool:
    return (path / "calculations").is_dir()


def read_setup(run: Path) -> SetupInfo:
    """Reads parameters.txt (info/setup/, or info/ for old results)."""
    info = SetupInfo()
    for candidate in (run / "info" / "setup" / "parameters.txt", run / "info" / "parameters.txt"):
        if candidate.exists():
            text = candidate.read_text(encoding="utf-8", errors="replace")
            info.found = True
            info.has_run_date = re.search(r"^Run date:", text, re.MULTILINE) is not None
            # "Computation time: 0h 00m 01s (0.9 s)"
            match = re.search(r"^Computation time:.*\(([0-9.]+) s\)", text, re.MULTILINE)
            if match:
                info.computation_seconds = float(match.group(1))
            match = re.search(r"^Git commit:\s*(\S+)", text, re.MULTILINE)
            if match:
                info.commit = match.group(1)
            break
    return info


def flatten(config: dict, prefix: str = "") -> Dict[str, object]:
    """{"solver": {"sizeX": 100}} -> {"solver.sizeX": 100}"""
    flat: Dict[str, object] = {}
    for key, value in config.items():
        if isinstance(value, dict):
            flat.update(flatten(value, prefix + key + "."))
        else:
            flat[prefix + key] = value
    return flat


def compare_pair(name: str, path_a: Path, path_b: Path, tolerance: float) -> PairResult:
    result = PairResult(name, path_a, path_b)
    for label, path in (("A", path_a), ("B", path_b)):
        if not is_run(path):
            result.notes.append(f"run {label} ({path}) does not exist or has no calculations/ folder")
    if result.notes:
        return result

    result.setup_a, result.setup_b = read_setup(path_a), read_setup(path_b)
    try:
        data_a, data_b = Dataset(path_a), Dataset(path_b)
    except (FileNotFoundError, ValueError) as error:
        result.notes.append(str(error))
        return result

    # A run that recorded its start but not its end crashed, was killed or is still running.
    # (Runs from before the run date was recorded have neither line; nothing can be said about them.)
    unfinished = [label for label, setup in (("A", result.setup_a), ("B", result.setup_b))
                  if setup.has_run_date and setup.computation_seconds is None]
    if unfinished:
        result.notes.append(f"run {' and '.join(unfinished)} did not finish (no 'Computation time' in parameters.txt)")
        return result

    config_a = flatten(data_a.configuration)
    config_b = flatten(data_b.configuration)
    grid_keys = ["solver.sizeX", "solver.sizeY", "solver.domain.x_left", "solver.domain.x_right",
                 "solver.domain.y_left", "solver.domain.y_right"]
    if any(config_a.get(key) != config_b.get(key) for key in grid_keys):
        result.verdict = "INCOMPATIBLE"
        result.notes.append("the grids differ: " + ", ".join(
            f"{key} {config_a.get(key)} vs {config_b.get(key)}" for key in grid_keys if config_a.get(key) != config_b.get(key)))
        return result
    for key in sorted(set(config_a) | set(config_b)):
        if config_a.get(key) != config_b.get(key):
            result.notes.append(f"config differs: {key} = {config_a.get(key)} (A) vs {config_b.get(key)} (B)")

    frames_a, frames_b = data_a.get_max_frame_number() + 1, data_b.get_max_frame_number() + 1
    if frames_a != frames_b:
        result.notes.append(f"different number of frames: {frames_a} (A) vs {frames_b} (B); comparing the first {min(frames_a, frames_b)}")

    nan_in = set()
    conc_min, conc_max = np.inf, -np.inf
    for frame in range(min(frames_a, frames_b)):
        _, _, phase_a = data_a.load(frame, "phase")
        _, _, conc_a = data_a.load(frame, "concentration")
        _, _, phase_b = data_b.load(frame, "phase")
        _, _, conc_b = data_b.load(frame, "concentration")

        identical = (phase_a.dtype == phase_b.dtype
                     and phase_a.tobytes() == phase_b.tobytes() and conc_a.tobytes() == conc_b.tobytes())
        for label, phase, conc in (("A", phase_a, conc_a), ("B", phase_b, conc_b)):
            if np.isnan(phase).any() or np.isnan(conc).any():
                nan_in.add(label)
        conc_min = min(conc_min, np.nanmin(conc_a))
        conc_max = max(conc_max, np.nanmax(conc_a))

        # The differences in double, so that float results do not lose digits when subtracted.
        result.rows.append(FrameRow(
            frame=frame,
            identical=identical,
            max_dphase=float(np.max(np.abs(phase_a.astype(np.float64) - phase_b.astype(np.float64)))),
            max_dconc=float(np.max(np.abs(conc_a.astype(np.float64) - conc_b.astype(np.float64)))),
            alpha_cells_a=int(np.count_nonzero(phase_a > 0.5)),
            alpha_cells_b=int(np.count_nonzero(phase_b > 0.5)),
        ))

    # Scale of the concentration differences: its range in run A. Without zirconium the
    # concentration stays constant (range 0), then its magnitude is used instead.
    conc_range = conc_max - conc_min
    if conc_range > 0:
        result.conc_scale = float(conc_range)
    elif np.isfinite(conc_max) and abs(conc_max) > 0:
        result.conc_scale = float(abs(conc_max))

    # A test is only meaningful when the solution changes during the run.
    _, _, first = data_a.load(0, "phase")
    first = first.copy()
    _, _, last = data_a.load(frames_a - 1, "phase")
    if np.array_equal(first, last):
        result.notes.append("the phase in run A does not change between the first and the last frame, so this comparison says little")

    if nan_in:
        result.notes.append(f"NaN in run {' and '.join(sorted(nan_in))}")
        result.verdict = "FAILED"
    elif frames_a == frames_b and all(row.identical for row in result.rows):
        result.verdict = "IDENTICAL"
    elif (frames_a == frames_b
          and result.max_dphase() <= tolerance
          and result.max_dconc_relative() <= tolerance
          and all(row.alpha_cells_a == row.alpha_cells_b for row in result.rows)):
        result.verdict = "CLOSE"
    else:
        result.verdict = "DIFFERENT"
    return result


def number(value: Optional[float]) -> str:
    if value is None:
        return "-"
    return "0" if value == 0 else f"{value:.2e}"


def seconds(setup: SetupInfo) -> str:
    return "-" if setup.computation_seconds is None else f"{setup.computation_seconds:.1f}"


def print_pair(result: PairResult) -> None:
    print(f"=== {result.name} ===")
    for label, path, setup in (("A", result.path_a, result.setup_a), ("B", result.path_b, result.setup_b)):
        details = []
        if setup.commit:
            details.append(f"commit {setup.commit}")
        if setup.computation_seconds is not None:
            details.append(f"computed in {setup.computation_seconds:.1f} s")
        print(f"{label}: {path}" + (f"   ({', '.join(details)})" if details else ""))

    if result.rows and result.verdict != "IDENTICAL":
        print(f"{'frame':>6}  {'identical':>9}  {'max|dphase|':>11}  {'max|dconc|/range':>16}  {'alpha cells A -> B':>18}")
        for row in result.rows:
            cells = f"{row.alpha_cells_a} -> {row.alpha_cells_b}"
            print(f"{row.frame:>6}  {'yes' if row.identical else 'no':>9}  {number(row.max_dphase):>11}  "
                  f"{number(row.max_dconc / result.conc_scale):>16}  {cells:>18}")
    elif result.rows:
        print(f"All {len(result.rows)} frames are bitwise identical.")
    print(f"Verdict: {result.verdict}")
    for note in result.notes:
        print(f"  - {note}")
    print()


def print_summary(results: List[PairResult], tolerance: float) -> None:
    width = max(len("case"), *(len(r.name) for r in results))
    print(f"Summary (tolerance {tolerance:g}):")
    print(f"{'case':<{width}}  {'verdict':<12}  {'max|dphase|':>11}  {'max|dconc|/range':>16}  {'time A [s]':>10}  {'time B [s]':>10}")
    for r in results:
        print(f"{r.name:<{width}}  {r.verdict:<12}  {number(r.max_dphase()):>11}  {number(r.max_dconc_relative()):>16}  "
              f"{seconds(r.setup_a):>10}  {seconds(r.setup_b):>10}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Compares the results of two runs (or two folders of runs) frame by frame.")
    parser.add_argument("run_a", type=Path, help="first run folder, or a folder of runs (the baseline)")
    parser.add_argument("run_b", type=Path, help="second run folder, or a folder of runs")
    parser.add_argument("--tol", type=float, default=DEFAULT_TOLERANCE,
                        help=f"largest difference for the verdict CLOSE (default {DEFAULT_TOLERANCE:g}); "
                             "absolute for the phase, relative to the concentration range for the concentration")
    args = parser.parse_args()

    for path in (args.run_a, args.run_b):
        if not path.is_dir():
            print(f"Folder {path} does not exist.", file=sys.stderr)
            return 2

    if is_run(args.run_a) or is_run(args.run_b):
        results = [compare_pair(f"{args.run_a.name} vs {args.run_b.name}", args.run_a, args.run_b, args.tol)]
    else:
        names = sorted({p.name for folder in (args.run_a, args.run_b) for p in folder.iterdir() if is_run(p)})
        if not names:
            print(f"Neither {args.run_a} nor {args.run_b} is a run folder or contains run folders.", file=sys.stderr)
            return 2
        results = [compare_pair(name, args.run_a / name, args.run_b / name, args.tol) for name in names]

    for result in results:
        print_pair(result)
    if len(results) > 1:
        print_summary(results, args.tol)
        print()
    counts = ", ".join(f"{sum(r.verdict == v for r in results)} {v}" for v in VERDICTS)
    print(f"Overall: {counts}")
    return 0 if all(r.verdict in ("IDENTICAL", "CLOSE") for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
