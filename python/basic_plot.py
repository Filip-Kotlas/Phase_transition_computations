"""Plots elementary plots of the calculations: selected frames and animations."""
import argparse
from cli import (add_common_arguments, frame_number, use_headless_backend, load_dataset,
                 decide_drawing, build_plotters, resolve_frames, render)
from plotter import save_animations

DEFAULT_FRAMES = [0, 1, 5, 10, 25, 50, 75, 100, 250, 500, 1000, "last"]


def parse_arguments(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plots selected frames and creates animations of the "
                                                 "computed results both in 2D and 3D.")
    add_common_arguments(parser)
    parser.add_argument("--format",
                        type=str,
                        choices=["auto", "mp4", "gif"],
                        help="Format of the animations. auto = mp4 if ffmpeg is available, otherwise gif.",
                        default="auto")
    parser.add_argument("--frames",
                        type=frame_number,
                        nargs="+",
                        metavar="FRAME",
                        help="Frames to be saved as images, numbers or 'last' for the last frame "
                             "(default: " + " ".join(str(frame) for frame in DEFAULT_FRAMES) + "). "
                             "Frames that do not exist are skipped.",
                        default=None)
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = parse_arguments(argv)
    use_headless_backend()

    dataset = load_dataset(args.name)
    plan = decide_drawing(dataset.configuration, args.type)
    plotters = build_plotters(dataset, plan)
    if not plotters:
        print("Nothing to plot.")
        return

    frames = resolve_frames(DEFAULT_FRAMES if args.frames is None else args.frames,
                            dataset.get_max_frame_number(),
                            report_unavailable=args.frames is not None)
    for frame in frames:
        for plotter in plotters:
            render(plotter, frame, save=True)

    save_animations(plotters, args.format)


if __name__ == "__main__":
    main()
