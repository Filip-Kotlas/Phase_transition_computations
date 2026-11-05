"""Shared parts of the command line scripts basic_plot.py and plot.py.

Both scripts have to be run from the root of the repository, because the results are looked up
in ./results/<name>.
"""
from dataclasses import dataclass
from pathlib import Path
from typing import List, Union
import argparse
import matplotlib
from dataset import Dataset
from plotter import Plotter, BoundaryPlotter2D, SurfacePlotter, CutPlotter

RESULTS_FOLDER = Path("results")
PLOT_TYPES = ["phase", "conc", "both"]
# Force terms for which the concentration is not computed (its RHS is zero), so there is nothing to draw.
FORCE_TERMS_WITHOUT_CONCENTRATION = ("radial", "constant")
# Initial conditions of the radial force term with a known analytic solution.
ANALYTIC_INITIAL_CONDITIONS = ("WulffShape", "ConstantCircle", "Circle")


def use_headless_backend() -> None:
    """Switches to the backend that only renders to files (no window, works without a display)."""
    matplotlib.use("Agg")


def frame_number(text: str) -> Union[int, str]:
    """Type of the command line arguments with a frame: a non-negative integer or the word "last"."""
    if text == "last":
        return "last"
    try:
        number = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError(f"'{text}' is not a frame number or 'last'.") from None
    if number < 0:
        raise argparse.ArgumentTypeError(f"Frame number cannot be negative, got {number}.")
    return number


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Adds the arguments --name and --type, which both scripts have."""
    parser.add_argument("--name",
                        type=str,
                        help="Name of the folder with results. Has to be in the results folder.",
                        default="Results")
    parser.add_argument("--type",
                        type=str,
                        choices=PLOT_TYPES,
                        help="Choose what will be plotted: the phase, the concentration or both.",
                        default="both")


def load_dataset(name: str) -> Dataset:
    """Opens the results folder with the given name."""
    results_path = RESULTS_FOLDER / name
    if not results_path.exists():
        raise FileNotFoundError(f"Složka {results_path} neexistuje!")
    return Dataset(results_path)


@dataclass
class DrawingPlan:
    """What is going to be drawn."""
    draw_phase: bool
    draw_concentration: bool
    draw_analytic_boundary: bool


def decide_drawing(configuration: dict, plot_type: str) -> DrawingPlan:
    """
    Decides from the requested type and the configuration of the run what is going to be drawn.

    The concentration is drawn only if it was computed (it is not for the radial and constant force
    terms). The analytic boundary is drawn only for the cases with a known analytic solution.
    """
    force_type = configuration["force_term"]["type"]
    initial_condition = configuration["init_cond"]["initial_condition"]

    draw_concentration = plot_type in ("conc", "both")
    if draw_concentration and force_type in FORCE_TERMS_WITHOUT_CONCENTRATION:
        if plot_type == "conc":
            print(f"The concentration is not computed for the force term '{force_type}', "
                  "so there is nothing to draw.")
        draw_concentration = False

    draw_analytic_boundary = ((force_type == "constant" and configuration["force_term"]["size"] == 0)
                              or (force_type == "radial"
                                  and initial_condition in ANALYTIC_INITIAL_CONDITIONS))

    return DrawingPlan(draw_phase=plot_type in ("phase", "both"),
                       draw_concentration=draw_concentration,
                       draw_analytic_boundary=draw_analytic_boundary)


def build_plotters(dataset: Dataset, plan: DrawingPlan) -> List[Plotter]:
    """Creates the plotters needed by the plan (2D plot, 3D surface and cut for the phase and the concentration)."""
    plotters: List[Plotter] = []
    if plan.draw_phase:
        plotters += [BoundaryPlotter2D(dataset, False, True, plan.draw_analytic_boundary, "phase"),
                     SurfacePlotter(dataset, "phase"),
                     CutPlotter(dataset, "phase", "x", 0.5, False)]
    if plan.draw_concentration:
        plotters += [BoundaryPlotter2D(dataset, True, False, False, "concentration"),
                     SurfacePlotter(dataset, "concentration"),
                     CutPlotter(dataset, "concentration", "x", 0.5, False)]
    return plotters


def resolve_frames(frames: List[Union[int, str]], max_frame: int, report_unavailable: bool) -> List[int]:
    """
    Turns the requested frames into the list of the available frame numbers.

    "last" becomes the number of the last frame, duplicates are dropped (the order is kept) and the
    frames that do not exist are skipped.

    Args:
        frames: Requested frames, numbers or "last".
        max_frame: Number of the last available frame.
        report_unavailable: If True, the skipped frames are reported.
    """
    resolved: List[int] = []
    unavailable: List[int] = []
    for frame in frames:
        number = max_frame if frame == "last" else frame
        if number > max_frame:
            unavailable.append(number)
        elif number not in resolved:
            resolved.append(number)
    if unavailable and report_unavailable:
        print(f"Skipping unavailable frames {unavailable} (the last frame is {max_frame}).")
    return resolved


def render(plotter: Plotter, frame: int, save: bool) -> None:
    """Saves the frame to a file, or shows it in a window."""
    if save:
        plotter.save_frame(frame)
    else:
        plotter.show_frame(frame)
