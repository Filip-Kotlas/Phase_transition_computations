from plotter import ScatterPlotter2D, SurfacePlotter

input_folder = "Results_rough_mesh_a_10"
dimension_2 = False

plotter = None
dimension_3 = not dimension_2
if dimension_2:
    plotter = ScatterPlotter2D(input_folder, False, True, True)
if dimension_3:
    plotter = SurfacePlotter(input_folder, True, False, False)

plotter.save_animation()
