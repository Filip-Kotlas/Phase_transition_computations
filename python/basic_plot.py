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

scatter_plotter = ScatterPlotter2D(results_path, False, True, False)
surface_plotter = SurfacePlotter(results_path, True, False, False)

scatter_plotter.save_frame(0)
scatter_plotter.save_frame(10)
scatter_plotter.save_frame(100)
scatter_plotter.save_frame(300)

surface_plotter.save_frame(0)
surface_plotter.save_frame(10)
surface_plotter.save_frame(100)
surface_plotter.save_frame(300)

scatter_plotter.save_animation()
surface_plotter.save_animation()