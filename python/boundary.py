"""The boundary of the phase (the level 0.5): the numerical one, the analytic one and their saving.

Nothing here draws. The plotters take the curves from the class Boundary and the function
write_boundary_files saves them for all the frames of a computation.
"""
from contextlib import ExitStack
from pathlib import Path
from typing import List, Tuple
import math
import numpy as np
import contourpy
from dataset import Dataset

NUMERICAL_FILE = "interface_num.crv"
ANALYTIC_FILE = "interface_anal.crv"

# The initial conditions whose boundary has an analytic solution: the circle and the Wulff shape. The circle
# is either sharp (ConstantCircle) or has a linear interface (Circle), in both cases the phase is 0.5
# exactly at the radius, so the boundary is the same.
CIRCLE_INITIAL_CONDITIONS = ("ConstantCircle", "Circle")
ANALYTIC_INITIAL_CONDITIONS = CIRCLE_INITIAL_CONDITIONS + ("WulffShape",)


def join_curves(curves: List[np.ndarray], separate_by_nan: bool=False) -> Tuple[np.ndarray, np.ndarray]:
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


class Boundary():
    """The boundary of the phase of one computation (a Dataset), as lists of curves.

    Both kinds of the curves are lists of arrays of the shape (number of points, 2), the columns are x
    and y. A closed curve ends with its first point again. A curve that reaches the border of
    the domain is open.
    """
    def __init__(self, dataset: Dataset) -> None:
        self.dataset = dataset
        self.configuration = dataset.configuration

        # In the data x changes fastest, so the shape of the grid is (sizeY, sizeX): the first index
        # is the position in y, the second one the position in x.
        self.grid_shape = (self.configuration["solver"]["sizeY"],
                           self.configuration["solver"]["sizeX"])
        x, y, _ = dataset.load(0, "phase")
        # Coordinates of the grid lines. They are copied, so the whole frame is not kept in memory by them.
        self.grid_x_lines = np.array(x.reshape(self.grid_shape)[0])
        self.grid_y_lines = np.array(y.reshape(self.grid_shape)[:, 0])

    def numerical_curves(self, frame: int, threshold: float=0.5) -> List[np.ndarray]:
        """
        The boundary of the computed phase in the given frame, i.e. its level "threshold".

        The phase is interpolated linearly along the edges of the grid. The points of the curves are the points
        on the edges where the interpolated phase equals the threshold. The curves are found by contourpy,
        which returns the points already ordered along every curve.
        """
        _, _, phase = self.dataset.load(frame, "phase")
        generator = contourpy.contour_generator(self.grid_x_lines, self.grid_y_lines,
                                                phase.reshape(self.grid_shape))
        return generator.lines(threshold)

    def analytic_curves(self, frame: int) -> List[np.ndarray]:
        """
        The boundary of the analytic solution in the given frame, as a list of curves like numerical_curves.

        The points of the curves are the exact intersections of the analytic boundary with the lines of the grid
        (x = const and y = const through the grid points). So they lie on the same edges of the grid as the points
        of the numerical boundary and they are spaced in the same way. Only the part of the boundary inside
        the domain is returned. The list is empty if there is no shape.
        """
        t = frame * self.dataset.time_step
        smallest_step = min(np.diff(self.grid_x_lines).min(), np.diff(self.grid_y_lines).min())

        # The exact boundary is replaced by a polyline whose segments are much shorter than the grid step,
        # so the error of this replacement is negligible.
        number_of_segments = 2000
        while True:
            curve_x, curve_y = analytic_polyline(self.configuration, t, number_of_segments)
            if len(curve_x) == 0:
                return []
            if np.hypot(np.diff(curve_x), np.diff(curve_y)).max() <= smallest_step / 50:
                break
            number_of_segments *= 2
            if number_of_segments > 100_000_000:
                raise RuntimeError("The analytic boundary cannot be approximated by a fine enough polyline.")
        return crossings_with_grid_lines(curve_x, curve_y, self.grid_x_lines, self.grid_y_lines)


def analytic_polyline(configuration: dict, t: float, number_of_segments: int=1000) -> Tuple[np.ndarray, np.ndarray]:
    """
    The analytic boundary at the time t as a closed polyline with number_of_segments segments (x and y coordinates).

    The shape is the circle (initial condition ConstantCircle or Circle) or the Wulff shape (WulffShape).
    Its squared radius changes linearly in time, it is empty if the radius is gone or if the force term
    has no analytic solution.
    """
    initial_condition = configuration["init_cond"]["initial_condition"]
    if initial_condition not in ANALYTIC_INITIAL_CONDITIONS:
        raise ValueError("Wrong initial condition for analytic solution.")

    # The radius is relative to the shorter side of the domain, the same as in the solver (Parameters.cu).
    domain = configuration["solver"]["domain"]
    smaller_side = min(domain["x_right"] - domain["x_left"],
                       domain["y_right"] - domain["y_left"])
    r0 = configuration["init_cond"]["radius_proportion"] * smaller_side

    force_term = configuration["force_term"]
    if force_term["type"] == "constant" and force_term["size"] == 0:
        r_squared = r0**2 - 2*t
    elif force_term["type"] == "radial":
        r_squared = r0**2 + 2*t
    else:
        return np.array([]), np.array([])

    # There is no shape.
    if r_squared < 0:
        return np.array([]), np.array([])
    r = math.sqrt(r_squared)

    offset_x = (domain["x_left"] + domain["x_right"])/2
    offset_y = (domain["y_left"] + domain["y_right"])/2
    theta = np.linspace(0, 2*math.pi, num=number_of_segments, endpoint=False)
    if initial_condition in CIRCLE_INITIAL_CONDITIONS:
        boundary_x = r * np.cos(theta) + offset_x
        boundary_y = r * np.sin(theta) + offset_y
    else:
        # The Wulff shape with the anisotropy psi(theta) = 1 + A cos(m (theta - theta_0)).
        anisotropy = configuration["anisotropy"]
        theta_0, m, A = anisotropy["theta_0"], anisotropy["m"], anisotropy["A"]
        psi = 1 + A*np.cos(m*(theta - theta_0))
        der_psi = -m*A*np.sin(m*(theta - theta_0))
        boundary_x = r*(psi*np.cos(theta) - der_psi*np.sin(theta)) + offset_x
        boundary_y = r*(psi*np.sin(theta) + der_psi*np.cos(theta)) + offset_y

    # The curve is closed by its first point.
    return np.append(boundary_x, boundary_x[0]), np.append(boundary_y, boundary_y[0])


def crossings_with_grid_lines(curve_x, curve_y, grid_x, grid_y) -> List[np.ndarray]:
    """
    Finds the points where a closed polyline crosses the lines of the grid, ordered along the polyline.

    The polyline has to be closed (its last point is its first point) and its segments have to be shorter
    than the grid step. Only the points inside the domain are kept. If the whole polyline is inside
    the domain, the result is one closed curve (it ends with its first point). Otherwise the result are open
    curves, one for every part of the polyline that lies inside the domain.

    Args:
        curve_x, curve_y: Coordinates of the points of the polyline.
        grid_x, grid_y: Coordinates of the vertical (x = const) and horizontal (y = const) lines of the grid.
    """
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


def write_boundary_files(dataset: Dataset, include_analytic: bool) -> None:
    """
    Saves the boundary of all the frames of the computation into the folder "info" of its results.

    The numerical boundary goes to interface_num.crv, the analytic one (if include_analytic) to
    interface_anal.crv. The files are always written from scratch, the frames follow one after another.
    The file has one block for every frame: the line "# Time level T = ...", the points "x y" and the line "##".
    """
    boundary = Boundary(dataset)
    info_folder = Path(dataset.folder) / "info"
    info_folder.mkdir(parents=True, exist_ok=True)

    outputs = [(info_folder / NUMERICAL_FILE, boundary.numerical_curves)]
    if include_analytic:
        outputs.append((info_folder / ANALYTIC_FILE, boundary.analytic_curves))

    with ExitStack() as stack:
        files = [stack.enter_context(path.open("w", encoding="utf-8")) for path, _ in outputs]
        for frame in range(dataset.get_max_frame_number() + 1):
            for file, (_, curves_of_frame) in zip(files, outputs):
                boundary_x, boundary_y = join_curves(curves_of_frame(frame))
                file.write(f"# Time level T =  {dataset.time(frame):.6f}\n")
                for x, y in zip(boundary_x, boundary_y):
                    file.write(f"    {x:.17g}   {y:.17g}\n")
                file.write("##\n\n")

    for path, _ in outputs:
        print("Boundary saved at:", path)
