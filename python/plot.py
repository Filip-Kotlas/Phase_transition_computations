from pathlib import Path
import argparse
from plotter import ScatterPlotter2D, SurfacePlotter

parser = argparse.ArgumentParser(description="Plots a given frame of the computed results both in \
                                 2D and 3D.")
parser.add_argument("name",
                    type=str,
                    help="Name of the folder with results. Has to be in results folder.",
                    default="Results")
parser.add_argument("frame",
                    type=int,
                    help="Number of the frame to be shown.",
                    default=0)
args = parser.parse_args()


results_path = Path()
results_path = results_path.parent / "results" / args.name

#plotter = ScatterPlotter2D(results_path, False, True, True)
#plotter.show_frame(args.frame)

plotter = SurfacePlotter(results_path, True, False, False)
plotter.show_frame(args.frame)
