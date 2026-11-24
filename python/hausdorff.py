"""Hausdorff distance of the interface (the level 0.5 of the phase) of a run from a reference run or
from the analytic solution, frame by frame.

Usage:
    python3 python/hausdorff.py RUN REFERENCE [--spacing S] [--output FILE]
    python3 python/hausdorff.py RUN --analytic [--spacing S] [--output FILE]

The Python counterpart of hausdorf.pas from Measurement/ (M. Benes), with shapely (GEOS) instead of
the Pascal sampling. In every frame:
    1. The interface of the run is the level 0.5 of its phase, found by contourpy (Boundary in
       boundary.py, the same curves as in the plots): one or more polylines.
    2. The interface it is compared with is either that of the reference run (found the same way)
       or, with --analytic, the exact interface of the verification cases of allen_cahn (the circle
       or the Wulff shape of analytic_polyline in boundary.py), replaced by a polygon whose distance
       from the exact curve is negligible (ANALYTIC_TOLERANCE) and cut off at the edges of the domain.
    3. The Hausdorff distance of the two sets of curves A and B is
           d_H = max( max over a in A of dist(a, B), max over b in B of dist(b, A) ).
       It is computed the way shapely.hausdorff_distance (the discrete Hausdorff distance of GEOS)
       computes it: the curves are sampled (shapely.segmentize splits every segment into equal parts
       at most delta = S * h apart, h = max(hx, hy) of the run) and the distance of every sample point
       of one curve from the whole other polyline (its segments, not only its vertices) is measured,
       in both directions. The result is the same as
           shapely.hausdorff_distance(shapely.segmentize(A, delta), shapely.segmentize(B, delta)),
       but the nearest segment is found through an STRtree, which is much faster (GEOS compares every
       point with every segment).
The sampling can only underestimate: the distance of a point from a curve changes at most as fast as
the point moves, so the true Hausdorff distance of the polylines is at most d_H + delta/2. To show
how far the result is from that, every frame is computed once more with the samples 2 * delta apart;
the report gives the largest change and warns when it is above CHECK_WARNING of the distance.

The frames are then combined into norms in time, as in lp_errors.py:
    L-infinity in time = max over the frames,
    L2 in time         = sqrt( integral of d_H(t)^2 dt ), the trapezoidal rule over the frames.
A frame where neither curve has an interface counts as distance 0; a frame where only one of them
has one is an error (the distance is not defined; hausdorf.pas wrote 10000 there).

Against a reference run, the checks of lp_errors.py apply (the same domain, times, number of frames
and model; other differences only warn). With --analytic, the run has to be an allen_cahn run with
an analytic solution: the initial condition ConstantCircle, Circle or WulffShape (not from a file)
and the force term constant with the size 0 (a shrinking shape, r^2 = r0^2 - 2t) or radial (a growing
one, r^2 = r0^2 + 2t), where t is the time since the initial condition.
Lines of the report that are not data rows start with "#". Exit code: 0 on success, 1 when the
distance cannot be computed, 2 for wrong arguments.
"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional
import argparse
import sys
import numpy as np
import shapely
from boundary import Boundary, analytic_polyline, ANALYTIC_INITIAL_CONDITIONS, CIRCLE_INITIAL_CONDITIONS
from dataset import Dataset, force_term, model_type
from lp_errors import RunError, cells_of, check_compatible, l2_in_time, max_in_time, open_run, steps_of

# The default distance of the sample points, as a fraction of h = max(hx, hy) of the run (--spacing).
DEFAULT_SPACING = 0.01
# The check (the samples 2 * delta apart) warns when it changes the distance by more than this
# fraction of the distance.
CHECK_WARNING = 0.01
# The polygon that replaces the analytic curve may differ from it by at most this fraction of the
# shorter side of the domain (estimated by the distance of the midpoints of the arcs from the chords).
ANALYTIC_TOLERANCE = 1e-10
MAX_ANALYTIC_SEGMENTS = 2**23


class Curves:
    """The interface in one frame: a list of polylines (arrays of the shape (number of points, 2))."""

    def __init__(self, polylines: List[np.ndarray]) -> None:
        # contourpy and the clipped analytic curve give polylines of at least two points.
        self.polylines = [np.asarray(polyline, dtype=np.float64) for polyline in polylines]
        self._tree = None

    def is_empty(self) -> bool:
        return len(self.polylines) == 0

    def samples(self, spacing: float) -> np.ndarray:
        """Points of the polylines at most spacing apart along every polyline: the vertices and the
        points that shapely.segmentize adds. A polyline whose segments are much shorter than spacing
        (the polygon of the analytic solution) keeps only every k-th vertex (and the last one), with
        k = floor(spacing / its longest segment), so the kept ones are still at most spacing apart."""
        thinned = []
        for polyline in self.polylines:
            k = int(spacing // np.max(np.hypot(*np.diff(polyline, axis=0).T)))
            if k >= 2:
                polyline = np.vstack([polyline[:-1:k], polyline[-1:]])
            thinned.append(polyline)
        lines = shapely.segmentize([shapely.linestrings(polyline) for polyline in thinned], spacing)
        return shapely.get_coordinates(lines)

    def tree(self) -> shapely.STRtree:
        """An STRtree of all the segments (built once)."""
        if self._tree is None:
            starts = np.vstack([polyline[:-1] for polyline in self.polylines])
            ends = np.vstack([polyline[1:] for polyline in self.polylines])
            self._tree = shapely.STRtree(shapely.linestrings(np.stack([starts, ends], axis=1)))
        return self._tree


def directed_distance(points: np.ndarray, curves: Curves) -> float:
    """max over the points of their distance from the curves (from the nearest segment)."""
    _, distances = curves.tree().query_nearest(shapely.points(points), return_distance=True)
    return float(np.max(distances))


def hausdorff_distance(a: Curves, b: Curves, spacing: float) -> float:
    """The Hausdorff distance of the curves, sampled at most spacing apart. Both must be non-empty."""
    return max(directed_distance(a.samples(spacing), b), directed_distance(b.samples(spacing), a))


@dataclass
class AnalyticSolution:
    """The analytic interface of an allen_cahn run as polygons, see analytic_curves()."""
    configuration: dict
    description: str
    max_segments: int = 0                 # the most segments of the polygon over the frames
    clipped_frames: List[int] = field(default_factory=list)   # frames where it leaves the domain

    def curves(self, frame: int, elapsed: float) -> Curves:
        """The analytic interface at the time elapsed since the initial condition."""
        domain = self.configuration["solver"]["domain"]
        smaller_side = min(domain["x_right"] - domain["x_left"], domain["y_right"] - domain["y_left"])
        tolerance = ANALYTIC_TOLERANCE * smaller_side
        # The polyline with 2n segments has the vertices of the one with n segments (the even ones,
        # the same angles) and the midpoints of its arcs (the odd ones): their distance from the chords
        # estimates how far the polygon with n segments is from the curve. The finer one is used.
        segments = 1000
        while True:
            x, y = analytic_polyline(self.configuration, elapsed, 2 * segments)
            if len(x) == 0:
                return Curves([])
            points = np.column_stack([x, y])
            even, odd = points[0::2], points[1::2]
            if point_segment_distances(odd, even[:-1], even[1:]).max() <= tolerance:
                break
            segments *= 2
            if 2 * segments > MAX_ANALYTIC_SEGMENTS:
                raise RunError(f"frame {frame}: the analytic interface cannot be replaced by a polygon "
                               f"with at most {MAX_ANALYTIC_SEGMENTS} segments")
        self.max_segments = max(self.max_segments, 2 * segments)

        inside = ((points[:, 0] >= domain["x_left"]) & (points[:, 0] <= domain["x_right"])
                  & (points[:, 1] >= domain["y_left"]) & (points[:, 1] <= domain["y_right"]))
        if np.all(inside):
            return Curves([points])
        # Only the part inside the domain, like the numerical interface, which ends at the edges.
        self.clipped_frames.append(frame)
        box = shapely.box(domain["x_left"], domain["y_left"], domain["x_right"], domain["y_right"])
        parts = shapely.get_parts(shapely.intersection(shapely.linestrings(points), box))
        # Points (the curve only touching an edge) are left out.
        return Curves([shapely.get_coordinates(part) for part in parts
                       if shapely.get_type_id(part) == shapely.GeometryType.LINESTRING])


def point_segment_distances(points: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """The distance of points[k] from the segment (starts[k], ends[k]), for every k."""
    direction = ends - starts
    length_squared = np.sum(direction**2, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        position = np.clip(np.sum((points - starts) * direction, axis=1) / length_squared, 0.0, 1.0)
    position = np.where(length_squared > 0, position, 0.0)
    return np.hypot(*(points - starts - position[:, np.newaxis] * direction).T)


def analytic_solution(run: Dataset) -> AnalyticSolution:
    """The analytic solution of the run; RunError if the run has none."""
    configuration = run.configuration
    model = model_type(configuration)
    init_cond = configuration["init_cond"]
    problems = []
    if model != "allen_cahn":
        problems.append(f"the model is {model}, not allen_cahn")
    else:
        force_type, force_size = force_term(configuration)
        if not (force_type == "radial" or (force_type == "constant" and force_size == 0)):
            problems.append(f"the force term is {force_type} with the size {force_size} "
                            "(only constant with the size 0 and radial have one)")
    if init_cond["initial_condition"] not in ANALYTIC_INITIAL_CONDITIONS:
        problems.append(f"the initial condition is {init_cond['initial_condition']} "
                        f"(only {', '.join(ANALYTIC_INITIAL_CONDITIONS)} have one)")
    if init_cond["init_cond_from_file"]:
        problems.append("the initial condition is read from a file")
    if problems:
        raise RunError(f"the run {run.folder} has no analytic solution: " + "; ".join(problems))
    shape = "circle" if init_cond["initial_condition"] in CIRCLE_INITIAL_CONDITIONS else "Wulff shape"
    force_type, _ = force_term(configuration)
    growth = "shrinking, r^2 = r0^2 - 2t" if force_type == "constant" else "growing, r^2 = r0^2 + 2t"
    return AnalyticSolution(configuration, f"the {shape} ({growth}, initial condition "
                                           f"{init_cond['initial_condition']}, force term {force_type})")


def numerical_curves(boundary: Boundary, frame: int, label: str) -> Curves:
    """The level 0.5 of the phase of the frame; RunError if the phase is not finite."""
    try:
        _, _, phase = boundary.dataset.load(frame, "phase")
    except ValueError as error:
        raise RunError(f"the {label} {boundary.dataset.folder}: {error}")
    if not np.all(np.isfinite(phase)):
        raise RunError(f"frame {frame} of the {label} {boundary.dataset.folder} contains NaN or inf")
    return Curves(boundary.numerical_curves(frame))


@dataclass
class HausdorffResult:
    run: Path
    reference: Optional[Path]           # None for the analytic solution
    reference_description: str
    model: str
    run_cells: tuple
    h: float                            # max(hx, hy) of the run
    spacing: float                      # delta, the largest distance of the sample points
    times: np.ndarray
    time_step: float
    distances: np.ndarray               # d_H in every frame (samples delta apart)
    check_distances: np.ndarray         # the same with the samples 2 * delta apart
    warnings: List[str] = field(default_factory=list)


def sampling_changes(distances: np.ndarray, check_distances: np.ndarray):
    """The change of the distance in every frame when the samples are 2 * delta apart: absolute and
    relative to the distance (0 where the distance is 0)."""
    changes = np.abs(distances - check_distances)
    with np.errstate(invalid="ignore", divide="ignore"):
        relative = np.where(distances > 0, changes / distances, 0.0)
    return changes, relative


def compute_distances(run_path: Path, reference_path: Optional[Path],
                      spacing_fraction: float = DEFAULT_SPACING) -> HausdorffResult:
    """The Hausdorff distance of the interfaces in every frame. reference_path None: from the
    analytic solution. Raises RunError."""
    if not spacing_fraction > 0:
        raise RunError("--spacing has to be positive")
    run = open_run(run_path, "run")
    analytic = None
    if reference_path is None:
        analytic = analytic_solution(run)
        warnings: List[str] = []
        reference_boundary = None
    else:
        reference = open_run(reference_path, "reference")
        warnings = check_compatible(run, reference)
        reference_boundary = Boundary(reference)
    run_boundary = Boundary(run)

    h = max(steps_of(run))
    spacing = spacing_fraction * h
    solver = run.configuration["solver"]
    frame_count = solver["frame_num"] + 1
    time_step = (solver["final_time"] - solver["initial_time"]) / solver["frame_num"]
    times = solver["initial_time"] + time_step * np.arange(frame_count)

    distances = np.zeros(frame_count)
    check_distances = np.zeros(frame_count)
    for frame in range(frame_count):
        run_curves = numerical_curves(run_boundary, frame, "run")
        if analytic is not None:
            other = analytic.curves(frame, time_step * frame)
            other_label = "the analytic solution"
        else:
            other = numerical_curves(reference_boundary, frame, "reference")
            other_label = "the reference"
        if run_curves.is_empty() and other.is_empty():
            continue                        # no interface in either: distance 0
        if run_curves.is_empty() or other.is_empty():
            which, missing = ("the run", other_label) if other.is_empty() else (other_label, "the run")
            raise RunError(f"frame {frame} (t = {times[frame]:.6g}): {which} has an interface, {missing} has none, "
                           "so the Hausdorff distance is not defined")
        distances[frame] = hausdorff_distance(run_curves, other, spacing)
        check_distances[frame] = hausdorff_distance(run_curves, other, 2 * spacing)

    if analytic is not None:
        description = (f"{analytic.description}, a polygon of up to {analytic.max_segments} segments "
                       f"(at most {ANALYTIC_TOLERANCE:.0e} x the shorter side of the domain from the curve)")
        if analytic.clipped_frames:
            warnings.append(f"the analytic interface leaves the domain in {len(analytic.clipped_frames)} frames "
                            f"(from frame {analytic.clipped_frames[0]}); only its part inside the domain is compared")
    else:
        ref_cells = cells_of(reference_boundary.dataset)
        description = f"{reference_path}   ({ref_cells[0]} x {ref_cells[1]} cells)"

    _, relative = sampling_changes(distances, check_distances)
    worst = int(np.argmax(relative))
    if relative[worst] > CHECK_WARNING:
        warnings.append(f"sampling at 2 x delta changes the distance by {100 * relative[worst]:.2g} % in frame {worst}; "
                        "a smaller --spacing may be needed")

    return HausdorffResult(run=run_path, reference=reference_path, reference_description=description,
                           model=model_type(run.configuration), run_cells=cells_of(run), h=h,
                           spacing=spacing, times=times, time_step=time_step, distances=distances,
                           check_distances=check_distances, warnings=warnings)


def format_report(result: HausdorffResult) -> str:
    """The report: a header, one row per frame and the norms in time. All lines but the rows start with #."""
    def number(value: float) -> str:
        return f"{value:.8e}"

    changes, relative = sampling_changes(result.distances, result.check_distances)
    worst = int(np.argmax(relative))
    lines = [f"# Run:       {result.run}   ({result.run_cells[0]} x {result.run_cells[1]} cells, h = {result.h:.6e})",
             f"# Reference: {result.reference_description}",
             f"# Model: {result.model}",
             f"# Sampling: at most delta = {result.spacing / result.h:g} h = {result.spacing:.3e} apart; the true distance "
             f"of the polylines is at most d_H + delta/2 = d_H + {result.spacing / 2:.3e}",
             f"# Check: with the samples 2 delta apart the distance changes by at most {number(changes.max())} "
             f"(relatively at most {100 * relative[worst]:.2g} %, frame {worst})"]
    lines += [f"# Warning: {warning}" for warning in result.warnings]
    lines.append(f"#{'frame':>5}{'time':>16}{'d_H':>16}{'d_H(2 delta)':>16}")
    for frame, time in enumerate(result.times):
        lines.append(f"{frame:>6}{number(time):>16}{number(result.distances[frame]):>16}"
                     f"{number(result.check_distances[frame]):>16}")
    lines += ["# Norms in time (L2 in time: trapezoidal rule over the frames):",
              f"#   L-infinity(0,T) = {number(max_in_time(result.distances))}",
              f"#   L2(0,T)         = {number(l2_in_time(result.distances, result.time_step))}"]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Hausdorff distance of the interface (phase = 0.5) of a run from a reference run or from "
                    "the analytic solution, frame by frame (the Python counterpart of hausdorf.pas).")
    parser.add_argument("run", type=Path, help="the run whose interface is measured")
    parser.add_argument("reference", type=Path, nargs="?", help="the reference run, usually a finer one")
    parser.add_argument("--analytic", action="store_true",
                        help="compare with the analytic solution instead of a reference run (allen_cahn only)")
    parser.add_argument("--spacing", type=float, default=DEFAULT_SPACING, metavar="S",
                        help=f"the largest distance of the sample points as a fraction of h = max(hx, hy) of the run "
                             f"(default {DEFAULT_SPACING})")
    parser.add_argument("--output", type=Path, help="also write the report into this file")
    args = parser.parse_args()
    if (args.reference is None) == (not args.analytic):
        parser.error("give either a reference run or --analytic")

    try:
        result = compute_distances(args.run, args.reference, args.spacing)
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
