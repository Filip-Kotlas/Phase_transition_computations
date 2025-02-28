from plotter import ScatterPlotter2D, SurfacePlotter

input_folders = ["Results_rough_mesh_a_10", "Results_finer_mesh_a_10"]
dimension_2 = True

"""
plotter = None
if dimension_2:
    plotter = ScatterPlotter2D(input_folder, False, True, True)
else:
    plotter = SurfacePlotter(input_folder, True, False, False)
"""

plotter = None
for folder in input_folders:

    plotter = ScatterPlotter2D(folder, False, True, True)
    plotter.save_frame(0)
    plotter.save_frame(10)
    plotter.save_frame(100)
    plotter.save_frame(300)
    plotter.save_animation()

    plotter = SurfacePlotter(folder, True, False, False)
    plotter.save_frame(0)
    plotter.save_frame(10)
    plotter.save_frame(100)
    plotter.save_frame(300)
    plotter.save_animation()
