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
                    help="Write what you want plotted. Options: conc, phase.",
                    default="phase")
args = parser.parse_args()


results_path = Path()
results_path = results_path.parent / "results" / args.name

if args.type == "phase":
    plotter = BoundaryPlotter2D(results_path, False, True, False, data_drawn="phase")
    plotter.show_frame(args.frame)
    plotter = SurfacePlotter(results_path, "phase")
    plotter.show_frame(args.frame)
elif args.type == "conc":
    plotter = CutPlotter(results_path, "concentration", "x", 0.5, False)
    plotter.show_frame(args.frame)
    plotter = SurfacePlotter(results_path, "concentration")
    plotter.show_frame(args.frame)
else:
    raise ValueError("Wrong argument given.")