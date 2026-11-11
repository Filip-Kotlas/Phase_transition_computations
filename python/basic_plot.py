"""Plots elementary plots of the calculations: selected frames and animations.

--name is a run folder in results/, or a folder that contains runs, e.g. a regression test
(results/MyTest with base/<case> and new/<case>). In the second case every run in it is plotted,
each into its own info/plots/.
"""
import argparse
import sys
from typing import List
from cli import (RESULTS_FOLDER, add_common_arguments, frame_number, use_headless_backend, load_dataset,
                 decide_drawing, build_plotters, resolve_frames, render)
from plotter import save_animations
from boundary import write_boundary_files

DEFAULT_FRAMES = [0, 1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, "last"]


def parse_arguments(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plots selected frames and creates animations of the "
                                                 "computed results both in 2D and 3D.")
    add_common_arguments(parser)
    parser.add_argument("--format",
                        type=str,
                        choices=["auto", "mp4", "gif"],
                        help="Format of the animations. auto = mp4 if ffmpeg is available, otherwise gif.",
                        default="auto")
    parser.add_argument("--frames",
                        type=frame_number,
                        nargs="+",
                        metavar="FRAME",
                        help="Frames to be saved as images, numbers or 'last' for the last frame "
                             "(default: " + " ".join(str(frame) for frame in DEFAULT_FRAMES) + "). "
                             "Frames that do not exist are skipped.",
                        default=None)
    parser.add_argument("--force",
                        action="store_true",
                        help="Plot even if some frames are missing (the existing animations, boundary files "
                             "and frame images of the same names are then overwritten).")
    return parser.parse_args(argv)


def find_runs(name: str) -> List[str]:
    """
    Names (relative to the results folder) of the runs to plot: the folder itself if it is a run
    (it has calculations/), otherwise all the runs anywhere inside it, sorted.
    A folder without any run is returned as it is, so that loading it reports the error.
    """
    folder = RESULTS_FOLDER / name
    if (folder / "calculations").is_dir() or not folder.is_dir():
        return [name]
    runs = sorted(path.parent for path in folder.rglob("calculations") if path.is_dir())
    if not runs:
        return [name]
    return [path.relative_to(RESULTS_FOLDER).as_posix() for path in runs]


def main(argv=None) -> None:
    args = parse_arguments(argv)
    use_headless_backend()

    runs = find_runs(args.name)
    if len(runs) == 1:
        plot_run(runs[0], args)
        return

    # A folder of runs: one failed run (e.g. a test case that crashed) should not stop the others.
    print(f"'{args.name}' contains {len(runs)} runs, plotting each of them.")
    failed = []
    for name in runs:
        print(f"\n=== {name} ===")
        try:
            plot_run(name, args)
        except Exception as error:   # report it and go on with the next run
            print(f"Plotting {name} failed: {type(error).__name__}: {error}")
            failed.append(name)
    if failed:
        print(f"\nPlotting failed for {len(failed)} of {len(runs)} runs: {', '.join(failed)}")
        sys.exit(1)


def plot_run(name: str, args: argparse.Namespace) -> None:
    """Plots the selected frames and the animations of one run."""
    dataset = load_dataset(name)
    if not dataset.is_complete() and not args.force:
        print(f"Only {len(dataset.file_paths)} of the {dataset.expected_frame_count()} frames of the run "
              f"'{name}' are present (deleted, or the run was interrupted).\n"
              "Plotting now would overwrite the existing animations, boundary files and frame images "
              "with incomplete ones, so nothing was done. Use --force to plot anyway.")
        return
    plan = decide_drawing(dataset.configuration, args.type)
    plotters = build_plotters(dataset, plan)
    if not plotters:
        print("Nothing to plot.")
        return

    if plan.draw_phase:
        write_boundary_files(dataset, plan.draw_analytic_boundary)

    frames = resolve_frames(DEFAULT_FRAMES if args.frames is None else args.frames,
                            dataset.get_max_frame_number(),
                            report_unavailable=args.frames is not None)
    for frame in frames:
        for plotter in plotters:
            render(plotter, frame, save=True)

    save_animations(plotters, args.format)


if __name__ == "__main__":
    main()
