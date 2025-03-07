from pathlib import Path
import argparse
from plotter import ScatterPlotter2D, SurfacePlotter

parser = argparse.ArgumentParser(description="Plots frames 0, 10, 100, and 300 and creates\
                                              animation of the computed results both in 2D and 3D.")
parser.add_argument("name",
                    type=str,
                    help="Name of the folder with results. Has to be in results folder.",
                    default="Results")
args = parser.parse_args()


results_path = Path()
results_path = results_path.parent / "results" / args.name

plotter = ScatterPlotter2D(results_path, False, True, True)
plotter.save_frame(0)
plotter.save_frame(10)
plotter.save_frame(100)
plotter.save_frame(300)
plotter.save_animation()

plotter = SurfacePlotter(results_path, True, False, False)
plotter.save_frame(0)
plotter.save_frame(10)
plotter.save_frame(100)
plotter.save_frame(300)
plotter.save_animation()
