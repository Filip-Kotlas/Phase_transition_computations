"""Plots one given frame of the computed results both in 2D and 3D."""
import argparse
from cli import (add_common_arguments, frame_number, use_headless_backend, load_dataset,
                 decide_drawing, build_plotters, resolve_frames, render)


def parse_arguments(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plots a given frame of the computed results both in "
                                                 "2D and 3D.")
    add_common_arguments(parser)
    parser.add_argument("--frame",
                        type=frame_number,
                        help="Number of the frame to be plotted, or 'last' for the last frame.",
                        default=0)
    parser.add_argument("--save",
                        action=argparse.BooleanOptionalAction,
                        help="Save the plotted frame as an image (default), or with --no-save show "
                             "it in a window.",
                        default=True)
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_arguments(argv)
    if args.save:
        use_headless_backend()

    dataset = load_dataset(args.name)
    frames = resolve_frames([args.frame], dataset.get_max_frame_number(), report_unavailable=True)
    if not frames:
        return

    plan = decide_drawing(dataset.configuration, args.type)
    for plotter in build_plotters(dataset, plan):
        render(plotter, frames[0], args.save)


if __name__ == "__main__":
    main()
