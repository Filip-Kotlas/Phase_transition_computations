"""Summarises the computation times of a benchmark (tests/benchmark.sh).

Usage:
    python3 python/compare_times.py BENCH_FOLDER

BENCH_FOLDER is the output folder of tests/benchmark.sh, e.g. results/MyBench. The script reads
BENCH_FOLDER/times.csv (one line per run) and prints for every case and version the median time
over the repetitions with the fastest and the slowest one, and the speedup

    speedup = median time of base / median time of new,

so 1.25x means that the new version needs 1/1.25 = 80 % of the time of the base version, and
0.80x that it is slower. A speedup marked "~" is within the noise: the ranges (fastest to slowest
run) of the two versions overlap.

The time is the solver's own computation time from parameters.txt (the time loop including frame
writing, without start-up, config and initial condition). Versions from before it was recorded
(commit 98e885c) do not write it; for a case where a run lacks it, the wall time of the whole
solver process is used for both versions instead, marked with "*".

The steps columns show the integration steps of the adaptive Kutta-Merson method ("Integration
steps" in parameters.txt): accepted ones, rejected attempts as "+M". A speedup can come from cheaper
steps or from fewer steps (e.g. when reordered arithmetic changes the error estimate), so the last
column is the speedup per accepted step:

    per-step speedup = (median time of base / steps of base) / (median time of new / steps of new).

Versions up to commit 2f0ab0f do not count the steps ("-"). The steps of the runs of one
version should be the same; if the repetitions differ, a note says so.

The last line is the geometric mean of the speedups (cluster.sh list shows it).
Exit code: 0 when every case has at least one successful run of both versions, 1 otherwise,
2 for wrong arguments.
"""
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Dict, List, Optional, Tuple
import argparse
import csv
import math
import re
import sys

VERSIONS = ["base", "new"]


@dataclass
class Run:
    repetition: int
    computation: Optional[float]   # None when the run did not record it
    wall: float
    exit_code: int


@dataclass
class Case:
    name: str
    runs: Dict[str, List[Run]] = field(default_factory=lambda: {v: [] for v in VERSIONS})

    def successful(self, version: str) -> List[Run]:
        return [run for run in self.runs[version] if run.exit_code == 0]

    def uses_wall_time(self) -> bool:
        return any(run.computation is None for v in VERSIONS for run in self.successful(v))

    def times(self, version: str) -> List[float]:
        wall = self.uses_wall_time()
        return [run.wall if wall else run.computation for run in self.successful(version)]

    def failed(self, version: str) -> int:
        return len(self.runs[version]) - len(self.successful(version))


def read_times(path: Path) -> List[Case]:
    cases: Dict[str, Case] = {}
    with path.open(newline="", encoding="utf-8") as file:
        for row in csv.DictReader(file):
            case = cases.setdefault(row["case"], Case(row["case"]))
            case.runs[row["version"]].append(Run(
                repetition=int(row["repetition"]),
                computation=float(row["computation_s"]) if row["computation_s"] else None,
                wall=float(row["wall_s"]),
                exit_code=int(row["exit_code"]),
            ))
    return list(cases.values())


def setup_line(folder: Path, key: str) -> Optional[str]:
    """The value of 'key: value' in parameters.txt of the run in folder, if there is one
    (in info/setup/, or in info/ for versions before commit 98e885c)."""
    parameters = folder / "info" / "setup" / "parameters.txt"
    if not parameters.exists():
        parameters = folder / "info" / "parameters.txt"
    if not parameters.exists():
        return None
    match = re.search(rf"^{key}:\s*(.+)$", parameters.read_text(encoding="utf-8", errors="replace"), re.MULTILINE)
    return match.group(1).strip() if match else None


def run_folder(bench: Path, version: str, case: str, repetition: int) -> Path:
    """Where tests/benchmark.sh puts a run: the first repetition complete, the others under repeats/."""
    if repetition == 1:
        return bench / version / case
    return bench / "repeats" / version / f"{case}-{repetition}"


Steps = Tuple[int, int]   # (accepted, rejected)


def read_steps(folder: Path) -> Optional[Steps]:
    """'Integration steps: 12345 (rejected: 67)' from parameters.txt, None when the run did not record it."""
    value = setup_line(folder, "Integration steps")
    match = re.fullmatch(r"(\d+) \(rejected: (\d+)\)", value) if value else None
    return (int(match.group(1)), int(match.group(2))) if match else None


def steps_cell(steps: Optional[Steps]) -> str:
    """'12345 (+67)': accepted (+rejected)."""
    return "-" if steps is None else f"{steps[0]} (+{steps[1]})"


def time_cell(times: List[float]) -> str:
    """'29.3 (29.1-29.6)': median (fastest-slowest)."""
    if not times:
        return "-"
    if len(times) == 1:
        return f"{times[0]:.1f}"
    return f"{median(times):.1f} ({min(times):.1f}-{max(times):.1f})"


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarises the computation times of a benchmark (tests/benchmark.sh).")
    parser.add_argument("folder", type=Path, help="output folder of tests/benchmark.sh, e.g. results/MyBench")
    args = parser.parse_args()

    times_file = args.folder / "times.csv"
    if not times_file.exists():
        print(f"{times_file} does not exist; is {args.folder} a benchmark?", file=sys.stderr)
        return 2
    cases = read_times(times_file)
    if not cases:
        print(f"{times_file} contains no runs.", file=sys.stderr)
        return 2

    # Header: what was compared, and where. Both versions ran on the same GPU.
    first = cases[0].name
    commits = {v: setup_line(args.folder / v / first, "Git commit") or "unknown" for v in VERSIONS}
    repetitions = max(len(case.runs[v]) for case in cases for v in VERSIONS)
    print(f"Benchmark: new ({commits['new']}) against base ({commits['base']}), {repetitions} run(s) of each version per case")
    gpu = setup_line(args.folder / "new" / first, "GPU")
    computer = setup_line(args.folder / "new" / first, "Computer")
    if gpu or computer:
        print(f"GPU: {gpu or 'unknown'}" + (f" on {computer}" if computer else ""))
    print()

    width = max(len("case"), *(len(case.name) + 1 for case in cases))   # + 1 for the marker *
    print(f"{'case':<{width}}  {'base [s]':<20}  {'new [s]':<20}  {'speedup':<9}  "
          f"{'steps base':<16}  {'steps new':<16}  per-step speedup")
    speedups: List[float] = []
    incomplete: List[str] = []
    notes: List[str] = []
    for case in cases:
        base, new = case.times("base"), case.times("new")
        marker = "*" if case.uses_wall_time() else ""
        if base and new:
            speedup = median(base) / median(new)
            speedups.append(speedup)
            noise = " ~" if max(min(base), min(new)) <= min(max(base), max(new)) else ""
            result = f"{speedup:.2f}x{noise}"
        else:
            incomplete.append(case.name)
            result = "-"

        # Steps of the first successful run of each version (they should not differ between repetitions).
        steps: Dict[str, Optional[Steps]] = {}
        for v in VERSIONS:
            counts = [read_steps(run_folder(args.folder, v, case.name, run.repetition)) for run in case.successful(v)]
            steps[v] = counts[0] if counts else None
            if len(set(counts)) > 1:
                notes.append(f"{case.name}: the {v} runs took different numbers of steps: "
                             + ", ".join(steps_cell(c) for c in counts))
        per_step = "-"
        if base and new and steps["base"] and steps["new"]:
            per_step = f"{(median(base) / steps['base'][0]) / (median(new) / steps['new'][0]):.2f}x"

        print(f"{case.name + marker:<{width}}  {time_cell(base):<20}  {time_cell(new):<20}  {result:<9}  "
              f"{steps_cell(steps['base']):<16}  {steps_cell(steps['new']):<16}  {per_step}")
        for v in VERSIONS:
            if case.failed(v):
                notes.append(f"{case.name}: {case.failed(v)} of {len(case.runs[v])} {v} run(s) failed; see solver.log "
                             f"in {v}/{case.name} (run 1) or repeats/{v}/{case.name}-K (run K)")

    print()
    print("Times: median of the runs, fastest-slowest in brackets. speedup = base / new (above 1: new is faster);")
    print("~ = within the noise (the ranges overlap).")
    print("steps = accepted integration steps (+rejected attempts); per-step speedup = speedup of the time per accepted step.")
    if any(case.uses_wall_time() for case in cases):
        print("* = wall time of the whole solver process: a run of this case did not record its computation time.")
    for note in notes:
        print(f"  - {note}")

    # The results should agree, otherwise the speedup compares different computations.
    comparison = args.folder / "comparison.txt"
    if comparison.exists():
        lines = comparison.read_text(encoding="utf-8", errors="replace").splitlines()
        if lines:
            verdicts = re.sub(r"^Overall:\s*", "", lines[-1])
            print(f"Results of base and new (python/compare_runs.py, details in comparison.txt): {verdicts}")
    print()

    if incomplete:
        print(f"No times of both versions for: {', '.join(incomplete)}")
    if speedups:
        mean = math.exp(sum(math.log(s) for s in speedups) / len(speedups))
        print(f"Overall: speedup {mean:.2f}x (base / new, geometric mean of {len(speedups)} case(s))")
    else:
        print("Overall: no case has times of both versions")
    return 1 if incomplete else 0


if __name__ == "__main__":
    sys.exit(main())
