from pathlib import Path
import argparse
from plotter import BoundaryPlotter2D, SurfacePlotter, CutPlotter

parser = argparse.ArgumentParser(description="Plots frames 0, 10, 100, and 300 and creates\
                                              animation of the computed results both in 2D and 3D.")
parser.add_argument("name",
                    type=str,
                    help="Name of the folder with results. Has to be in results folder.",
                    default="Results")
args = parser.parse_args()


results_path = Path()
results_path = results_path.parent / "results" / args.name

boundary_plotter = BoundaryPlotter2D(results_path, False, True, True, "phase")
phase_surface_plotter = SurfacePlotter(results_path, "phase")
concentration_surface_plotter = SurfacePlotter(results_path, "concentration")
concentration_cut_plotter = CutPlotter(results_path, "concentration", "x", 0.5, True)


boundary_plotter.save_frame(0)
boundary_plotter.save_frame(10)
boundary_plotter.save_frame(100)
boundary_plotter.save_frame(300)

phase_surface_plotter.save_frame(0)
phase_surface_plotter.save_frame(10)
phase_surface_plotter.save_frame(100)
phase_surface_plotter.save_frame(300)


concentration_surface_plotter.save_frame(0)
concentration_surface_plotter.save_frame(10)
concentration_surface_plotter.save_frame(100)
concentration_surface_plotter.save_frame(300)

concentration_cut_plotter.save_frame(0)
concentration_cut_plotter.save_frame(10)
concentration_cut_plotter.save_frame(100)
concentration_cut_plotter.save_frame(300)

boundary_plotter.save_animation()
concentration_cut_plotter.save_animation()
concentration_surface_plotter.save_animation()
phase_surface_plotter.save_animation()
