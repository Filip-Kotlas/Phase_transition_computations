"""Plots elementar plots of the calculations."""
from pathlib import Path
import argparse
from plotter import BoundaryPlotter2D, SurfacePlotter, CutPlotter, save_animations
from dataset import Dataset

parser = argparse.ArgumentParser(description="Plots frames 0, 10, 50, 100 and creates\
                                              animation of the computed results both in 2D and 3D.")
parser.add_argument("--name",
                    type=str,
                    help="Name of the folder with results. Has to be in results folder.",
                    default="Results")
parser.add_argument("--type",
                    type=str,
                    help="Choose what will be plotted. Options: conc, phase, both",
                    default="both")
parser.add_argument("--format",
                    type=str,
                    choices=["auto", "mp4", "gif"],
                    help="Format of the animations. auto = mp4 if ffmpeg is available, otherwise gif.",
                    default="auto")
args = parser.parse_args()

results_path = Path()
results_path = results_path.parent / "results" / args.name

if not results_path.exists():
    raise FileNotFoundError(f"Soubor {results_path} neexistuje!")
dataset = Dataset(results_path)
configuration = dataset.configuration

draw_concentration = args.type in ("conc", "both")
if configuration["force_term"]["type"] in ("radial", "constant"):
    draw_concentration = False

draw_anal_bound = False
if ((configuration["force_term"]["type"] == "constant"
     and configuration["force_term"]["size"] == 0) 
    or (configuration["force_term"]["type"] == "radial"
        and (configuration["init_cond"]["initial_condition"] == "WulffShape"
             or configuration["init_cond"]["initial_condition"] == "ConstantCircle"
             or configuration["init_cond"]["initial_condition"] == "Circle"))):
    draw_anal_bound = True

boundary_plotter = BoundaryPlotter2D(dataset, False, True, draw_anal_bound, "phase")
concentration_2D_plotter = BoundaryPlotter2D(dataset, True, False, False, "concentration")
phase_surface_plotter = SurfacePlotter(dataset, "phase")
concentration_surface_plotter = SurfacePlotter(dataset, "concentration")
concentration_cut_plotter = CutPlotter(dataset, "concentration", "x", 0.5, False)
phase_cut_plotter = CutPlotter(dataset, "phase", "x", 0.5, False)
for i in [0, 1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, boundary_plotter.get_max_frame_number()]:
    if args.type in ("phase", "both"):
        boundary_plotter.save_frame(i)
        phase_cut_plotter.save_frame(i)
        phase_surface_plotter.save_frame(i)

    if draw_concentration:
        concentration_surface_plotter.save_frame(i)
        concentration_cut_plotter.save_frame(i)
        concentration_2D_plotter.save_frame(i)

animated_plotters = []
if draw_concentration:
    animated_plotters += [concentration_cut_plotter, concentration_surface_plotter, concentration_2D_plotter]

if args.type in ("phase", "both"):
    animated_plotters += [boundary_plotter, phase_surface_plotter, phase_cut_plotter]

save_animations(animated_plotters, args.format)

if args.type not in ("phase", "conc", "both"):
    raise ValueError("Wrong argument for type.")