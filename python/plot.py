from pathlib import Path
import argparse
from plotter import BoundaryPlotter2D, SurfacePlotter, CutPlotter

parser = argparse.ArgumentParser(description="Plots a given frame of the computed results both in \
                                 2D and 3D.")
parser.add_argument("--name",
                    type=str,
                    help="Name of the folder with results. Has to be in results folder.",
                    default="Results")
parser.add_argument("--frame",
                    type=int,
                    help="Number of the frame to be shown.",
                    default=0)
parser.add_argument("--type",
                    type=str,
                    help="Write what you want plotted. Options: conc, phase, both.",
                    default="both")
parser.add_argument("--save",
                    type=bool,
                    help="If True, saves the plotted frame as png file. If False, shows the frame in a window.",
                    default=True)
args = parser.parse_args()


results_path = Path()
results_path = results_path.parent / "results" / args.name

if args.type in ("phase", "both"):
    plotter = BoundaryPlotter2D(results_path, False, True, False, "phase")
    if args.save:
        plotter.save_frame(args.frame)
    else:
        plotter.show_frame(args.frame)
    plotter = SurfacePlotter(results_path, "phase")
    if args.save:
        plotter.save_frame(args.frame)
    else:
        plotter.show_frame(args.frame)
    plotter = CutPlotter(results_path, "phase", "x", 0.5, False)
    if args.save:
        plotter.save_frame(args.frame)
    else:
        plotter.show_frame(args.frame)
if args.type in ("conc", "both"):
    plotter = BoundaryPlotter2D(results_path, True, False, False, "concentration")
    if args.save:
        plotter.save_frame(args.frame)
    else:
        plotter.show_frame(args.frame)
    plotter = CutPlotter(results_path, "concentration", "x", 0.5, False)
    if args.save:
        plotter.save_frame(args.frame)
    else:
        plotter.show_frame(args.frame)
    plotter = SurfacePlotter(results_path, "concentration")
    if args.save:
        plotter.save_frame(args.frame)
    else:
        plotter.show_frame(args.frame)
else:
    raise ValueError("Wrong argument given.")