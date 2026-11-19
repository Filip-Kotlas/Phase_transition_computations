"""Experimental order of convergence (EOC) of several runs against a reference run.

Usage:
    python3 python/eoc.py RUN1 RUN2 ... --reference REF [--eval-cells NX NY] [--output FILE] [--tex FILE]

The Python version of eoc.cpp from Measurement/ (M. Benes), with the errors computed directly from
the result folders instead of from a max.txt filled in by hand. For every run:
    - the errors against the reference come from lp_errors.py (compute_errors): for every measured
      field (phase, and concentration if the model computes it) four norms in time and space:
          L-infinity(0,T; L2), L-infinity(0,T; L-infinity)   (the max over the frames, as in eoc.cpp)
          L2(0,T; L2),         L2(0,T; L-infinity)           (trapezoidal rule over the frames)
    - h = max(hx, hy) of the run (as lp_dist_new_new.pas writes into max.txt),
    - the number of accepted integration steps and the computation time from parameters.txt
      ("-" when the run did not record them),
    - DoF = the number of unknowns of the solver, 2 * sizeX * sizeY (phase and concentration).
The runs have to be given from the coarsest to the finest (h strictly decreasing). For every pair of
neighbouring runs and every error e:
    EOC_i = log(e_i / e_{i-1}) / log(h_i / h_{i-1})
as in eoc.cpp. The first run has no EOC ("-"), and neither has a run whose error or the previous
error is 0 (eoc.cpp wrote 0 or NaN there).

The report has, for every field, a table of the errors and a table of the EOCs. --output writes it
into a file too, --tex writes the same tables as LaTeX (no packages or macros needed).
The checks of lp_errors.py apply to every run against the reference (the same domain, times, number
of frames and model), so all the runs are compatible with each other too.
Exit code: 0 on success, 1 when the runs cannot be compared, 2 for wrong arguments.
"""
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple
import argparse
import math
import sys
from compare_runs import read_setup
from compare_times import read_steps
from lp_errors import FieldErrors, RunError, compute_errors, l2_in_time, max_in_time, open_run, steps_of


@dataclass
class Norm:
    label: str                                       # in the text report
    latex: str                                       # in the LaTeX tables
    value: Callable[[FieldErrors, float], float]     # (errors of one field, time between frames) -> norm


NORMS = [
    Norm("Linf(0,T;L2)", r"$L^\infty(0,T;L^2)$", lambda errors, dt: max_in_time(errors.l2)),
    Norm("Linf(0,T;Linf)", r"$L^\infty(0,T;L^\infty)$", lambda errors, dt: max_in_time(errors.linf)),
    Norm("L2(0,T;L2)", r"$L^2(0,T;L^2)$", lambda errors, dt: l2_in_time(errors.l2, dt)),
    Norm("L2(0,T;Linf)", r"$L^2(0,T;L^\infty)$", lambda errors, dt: l2_in_time(errors.linf, dt)),
]


@dataclass
class RunRow:
    """One run: one row of the tables."""
    path: Path
    h: float
    cells: Tuple[int, int]
    reference_cells: Tuple[int, int]
    eval_cells: Tuple[int, int]
    dof: int
    steps: Optional[int]                 # accepted integration steps
    seconds: Optional[float]             # computation time
    errors: Dict[str, List[float]]       # field -> the errors in the order of NORMS
    warnings: List[str]


def eoc_column(errors: Sequence[float], hs: Sequence[float]) -> List[Optional[float]]:
    """EOC_i = log(e_i / e_{i-1}) / log(h_i / h_{i-1}); None for the first row and where an error is 0."""
    column: List[Optional[float]] = [None]
    for i in range(1, len(errors)):
        if errors[i] > 0 and errors[i - 1] > 0:
            column.append(math.log(errors[i] / errors[i - 1]) / math.log(hs[i] / hs[i - 1]))
        else:
            column.append(None)
    return column


def measure(run_paths: Sequence[Path], reference: Path,
            eval_cells: Optional[Sequence[int]] = None) -> List[RunRow]:
    """The errors and the other columns of every run. Raises RunError."""
    if len(run_paths) < 2:
        raise RunError("the EOC needs at least two runs")
    # The order is checked before the (possibly long) computation of the errors.
    hs = [max(steps_of(open_run(path, "run"))) for path in run_paths]
    for i in range(1, len(hs)):
        if not hs[i] < hs[i - 1]:
            raise RunError(f"give the runs from the coarsest to the finest: h = {hs[i - 1]:.6g} ({run_paths[i - 1]}) "
                           f"is followed by h = {hs[i]:.6g} ({run_paths[i]})")
    rows = []
    for path in run_paths:
        print(f"Computing the errors of {path} ...", file=sys.stderr)
        result = compute_errors(path, reference, eval_cells)
        size_x, size_y = result.run_cells[0] + 1, result.run_cells[1] + 1
        steps = read_steps(path)
        rows.append(RunRow(
            path=path, h=result.h, cells=result.run_cells, reference_cells=result.reference_cells,
            eval_cells=result.eval_cells,
            dof=2 * size_x * size_y,
            steps=steps[0] if steps else None,
            seconds=read_setup(path).computation_seconds,
            errors={name: [norm.value(errors, result.time_step) for norm in NORMS]
                    for name, errors in result.fields.items()},
            warnings=result.warnings))
    return rows


def eocs(rows: List[RunRow], field: str) -> List[List[Optional[float]]]:
    """The EOCs of one field: one column per norm (in the order of NORMS), one value per run."""
    hs = [row.h for row in rows]
    return [eoc_column([row.errors[field][k] for row in rows], hs) for k in range(len(NORMS))]


def optional(value, text: Callable[[object], str]) -> str:
    return "-" if value is None else text(value)


def format_report(rows: List[RunRow], reference: Path) -> str:
    lines = [f"Reference: {reference}"]
    for row in rows:
        lines.append(f"Run:       {row.path}   ({row.cells[0]} x {row.cells[1]} cells, "
                     f"evaluation grid {row.eval_cells[0]} x {row.eval_cells[1]})")
    for row in rows:
        lines += [f"Warning ({row.path.name}): {warning}" for warning in row.warnings]
    name_width = max(len("run"), *(len(row.path.name) for row in rows))

    for field in rows[0].errors:
        lines += ["", f"Errors of the {field}:"]
        header = (f"{'run':<{name_width}}  {'h':>12}  {'steps':>10}  {'DoF':>10}  "
                  + "  ".join(f"{norm.label:>14}" for norm in NORMS) + f"  {'time [s]':>10}")
        lines += [header, "-" * len(header)]
        for row in rows:
            lines.append(f"{row.path.name:<{name_width}}  {row.h:>12.6g}  {optional(row.steps, str):>10}  {row.dof:>10}  "
                         + "  ".join(f"{error:>14.6e}" for error in row.errors[field])
                         + f"  {optional(row.seconds, lambda s: f'{s:.1f}'):>10}")

        lines += ["", f"EOC of the {field}:"]
        header = f"{'run':<{name_width}}  {'h':>12}  " + "  ".join(f"{norm.label:>14}" for norm in NORMS)
        lines += [header, "-" * len(header)]
        columns = eocs(rows, field)
        for i, row in enumerate(rows):
            lines.append(f"{row.path.name:<{name_width}}  {row.h:>12.6g}  "
                         + "  ".join(f"{optional(column[i], lambda v: f'{v:.3f}'):>14}" for column in columns))
    return "\n".join(lines)


def latex_number(value: float) -> str:
    """1.234e-05 -> $1.234\\cdot 10^{-5}$"""
    if value == 0:
        return "$0$"
    mantissa, exponent = f"{value:.3e}".split("e")
    return f"${mantissa}\\cdot 10^{{{int(exponent)}}}$"


def format_latex(rows: List[RunRow]) -> str:
    reference_cells = rows[0].reference_cells
    lines = []
    for field in rows[0].errors:
        lines += [r"\begin{table}[ht]", r"\centering",
                  r"\begin{tabular}{|l|r|r|" + "l|" * len(NORMS) + "r|}", r"\hline",
                  r"$h$ & steps & DoF & " + " & ".join(norm.latex for norm in NORMS) + r" & time [s] \\",
                  r"\hline"]
        for row in rows:
            lines.append(f"{row.h:.6g} & {optional(row.steps, str)} & {row.dof} & "
                         + " & ".join(latex_number(error) for error in row.errors[field])
                         + f" & {optional(row.seconds, lambda s: f'{s:.1f}')} \\\\")
        lines += [r"\hline", r"\end{tabular}",
                  f"\\caption{{Errors of the {field} against the reference run with "
                  f"${reference_cells[0]} \\times {reference_cells[1]}$ cells.}}",
                  r"\end{table}", ""]

        lines += [r"\begin{table}[ht]", r"\centering",
                  r"\begin{tabular}{|l|" + "r|" * len(NORMS) + "}", r"\hline",
                  r"$h$ & " + " & ".join(norm.latex for norm in NORMS) + r" \\",
                  r"\hline"]
        columns = eocs(rows, field)
        for i, row in enumerate(rows):
            lines.append(f"{row.h:.6g} & "
                         + " & ".join(optional(column[i], lambda v: f"{v:.3f}") for column in columns) + r" \\")
        lines += [r"\hline", r"\end{tabular}", f"\\caption{{EOC of the {field}.}}", r"\end{table}", ""]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Experimental order of convergence of several runs against a reference run "
                    "(the Python version of eoc.cpp, with the errors of lp_errors.py).")
    parser.add_argument("runs", type=Path, nargs="+", help="the runs, from the coarsest to the finest")
    parser.add_argument("--reference", type=Path, required=True, help="the reference run, finer than all the runs")
    parser.add_argument("--eval-cells", type=int, nargs=2, metavar=("NX", "NY"),
                        help="cells of the evaluation grid for every run (default: for each run the least "
                             "common multiple of its cell counts and those of the reference)")
    parser.add_argument("--output", type=Path, help="also write the report into this file")
    parser.add_argument("--tex", type=Path, help="write the tables as LaTeX into this file")
    args = parser.parse_args()

    try:
        rows = measure(args.runs, args.reference, args.eval_cells)
    except RunError as error:
        print(f"Error: {error}", file=sys.stderr)
        return 1
    report = format_report(rows, args.reference)
    print(report)
    if args.output is not None:
        args.output.write_text(report + "\n", encoding="utf-8")
    if args.tex is not None:
        args.tex.write_text(format_latex(rows) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
