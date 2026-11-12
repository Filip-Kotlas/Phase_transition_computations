# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

CUDA/C++17 phase-field solver for a 2D phase transition (anisotropic Allen–Cahn-type equation for the phase `p`, coupled with diffusion of the Nb concentration `c` in a Zr–Nb alloy). It is research code for a thesis at FNSPE CTU Prague. Python scripts post-process and plot the results. Comments and messages are in English.

## Build & run

The code needs `nvcc` and the header-only **TNL** library (Template Numerical Library), which is expected in `~/.local/include`. Builds only work on a Linux machine with an NVIDIA GPU. The Makefile compiles with `-arch=native` (override with `make CUDA_ARCH=...`) and keeps objects in `build/sm_<compute capability>/`, so build on the machine that will run the binary. Jobs for the GPU cluster are submitted from a local machine with `slurm/cluster.sh` (see `slurm/run.sbatch`); from Windows PowerShell use the launcher `.\slurm\cluster.cmd`, which runs the script with Git's bash and Windows OpenSSH (the key lives in the Windows ssh-agent).

```bash
make                                  # build ./build/Phase-field
make run FOLDER=MyRun                 # build, run with config/config.json -> results/MyRun, then plot
make run FOLDER=MyRun DELETE_CALCS=1  # same, but keep only the first and last Result-*.npy afterwards
make run_all                          # run every comp_configurations/*.json (each one is copied over config/config.json)
make clean
```

- The binary is called as `Phase-field [result_folder [config_file]]`; the config defaults to `config/config.json` (`make run CONFIG=...` overrides it). `slurm/cluster.sh submit` uploads the sources of each run (including `config/config.json`) to `jobs/<FOLDER>/sources.tar.gz` on the cluster. The job unpacks them into the node's `/tmp`, builds them from scratch and runs `make run` there, with `results` linked to the shared results folder. A queued job therefore always runs the code it was submitted with, and jobs share no sources or build.
- Run everything from the repo root. The binary resolves its default config through a relative path, and the Python scripts look in `./results/<name>`.
- Header dependencies are tracked: nvcc writes `build/sm_XX/*.d` (`-MMD -MP`) and the Makefile includes them, so changing an `include/*.hpp` recompiles the objects that use it. System headers (TNL in `/usr/include`) are not tracked. `make clean` is still needed after changing compiler flags in the Makefile or after updating TNL (not on the cluster, where every job builds from scratch). (Before this, stale `.o` files after header changes caused hard-to-find bugs.)
- `make run_all` overwrites the tracked `config/config.json`. `comp_configurations/` is gitignored, apart from `.gitkeep`.
- There is no linter. There is a regression test and a benchmark (see below).

## Regression test

It answers "did my change alter the results?", not "are the results correct?" (that is the job of a comparison with the analytic solutions).

```bash
.\slurm\cluster.cmd test MyTest                     # working tree (uncommitted changes) vs HEAD
.\slurm\cluster.cmd test MyTest --against HEAD~1    # vs any commit, e.g. after committing the changes
.\slurm\cluster.cmd test MyTest --against A --new B # two commits: B is tested instead of the working tree
.\slurm\cluster.cmd status                          # while it runs: which build or run it is at
.\slurm\cluster.cmd log                             # the result; also in results/MyTest/report.txt
```
- `cluster.sh test` uploads the working tree (like `submit`) plus `git archive <commit>` as the baseline, and refuses when both have the same code (`config/config.json` does not count). So test **before** committing, or use `--against`.
- `--new COMMIT` (test and bench) uploads that commit as `new.tar.gz` and tests it instead of the working tree. The working tree is still uploaded as the harness: `run.sbatch`, `tests/`, the configs and `python/compare_*.py` always come from it (`harness_dir` in `tests/common.sh`), so commits from before the test existed can be compared too, as long as they accept `Phase-field RESULT_FOLDER CONFIG` and today's configs.
- The job (`run.sbatch` with `MODE=test`, time limit 1 h) runs `tests/regression.sh`: it builds both versions, runs every `tests/configs/*.json` with both on the same GPU (so GPU-type differences cannot appear) into `results/MyTest/{base,new}/<case>`, then runs `python/compare_runs.py`. Nothing is stored in the repo; the baseline is recomputed every time. The configs, scripts and comparison come from the new version.
- Cases (5 frames, seconds each): `circle_constant` (isotropic, constant force; 120×100 on [0, 1.5] × [0, 1], so hx ≠ hy and a mistake in the IC coordinates shows), `wulff_radial` (anisotropy, Wulff IC, radial force), `zirconium` (concentration), `dirichlet_mixed` (Dirichlet and Neumann for both fields), `from_file` (IC from a float32 `.npy` file with `scaling`, generated by `tests/prepare_input.py`), these four 100×100 on the unit square; `random_rectangle` (120×80 on [-0.5, 1.5] × [0.2, 1]; random values at every point from `random.npy`, without scaling, zirconium, all Neumann). The other cases have a constant phase near the edges, so a wrong boundary index would not show in them; `random_rectangle` catches it (`prepare_input.py` writes random values for any `random*.npy`). Keep at least one case with hx ≠ hy, otherwise swapped x and y go unnoticed. Each case should run a different code path; add a config for a new feature.
- Verdicts: IDENTICAL (bitwise), CLOSE (phase difference and concentration difference relative to its range ≤ `--tol`, default 1e-4, and no cell crossed 0.5), DIFFERENT, FAILED (missing, unfinished or NaN), INCOMPATIBLE (different grids). Cleanup should give IDENTICAL; reordered arithmetic (performance work) at best CLOSE. For reference, float vs double on the 200×200 zirconium run differed by 3.5e-4 in the phase after 100 frames.
- The test draws no plots. After `cluster.cmd fetch MyTest`, `python3 python/basic_plot.py --name MyTest` plots every run in it (it does that for any `--name` folder that is not a run itself but contains runs); one case: `--name MyTest/new/zirconium`.
- `python3 python/compare_runs.py results/A results/B` compares any two runs (or two folders of runs) the same way, e.g. two `submit` runs on different GPUs.

## Benchmark

It answers "how much faster is my change?". Built like the regression test: working tree vs a baseline commit, both on the same GPU, nothing stored in the repo.

```bash
.\slurm\cluster.cmd bench MyBench                          # working tree vs HEAD, 3 runs of each version per case
.\slurm\cluster.cmd bench MyBench --against HEAD~1 --repeat 5 --gpu rtx3060
.\slurm\cluster.cmd bench MyBench --against 5e9bccb~1 --new 5e9bccb   # two commits, e.g. an old change
.\slurm\cluster.cmd bench MyBench --long                   # one 800x800 case, one ~10 min run of each version
.\slurm\cluster.cmd log                                    # the result; also in results/MyBench/report.txt
```
- The job (`run.sbatch` with `MODE=bench`, time limit 1 h + 15 min per run above 3) runs `tests/benchmark.sh`: it builds both versions and runs every `tests/benchmark_configs/*.json` N times with each, alternating the order (base new, new base, …), so a drift on the node hits both alike. `tests/common.sh` holds the build and run functions shared with `regression.sh`.
- Time = the solver's own `Computation time` from `parameters.txt` (time loop including frame writing, without start-up and IC). Baselines older than 98e885c do not write it; then the wall time of the process is used for that case (marked `*`). `python/compare_times.py` reports median (fastest–slowest) per version and the speedup base/new (> 1: new is faster; `~` = the ranges overlap, i.e. within the noise), and the geometric mean over the cases.
- Output in `results/MyBench/`: `times.csv` (every run), `report.txt`, `{base,new}/<case>` (first run, complete) with `comparison.txt` from `compare_runs.py` (a speedup only counts if the results still agree), `repeats/` (the other runs, without frames).
- Cases (about half a minute each on an RTX 2070 SUPER, scaled from Bench_noise, where the runs of one case differed by ≤ 1 %): `zirconium_200` (the production config, shortened: the time per frame is constant during a run), `zirconium_400` (the same at 400×400), `wulff_radial_300` (anisotropy and radial force, no concentration; the Wulff shape shrinks and would vanish at t ≈ 0.055, so its length comes from the resolution, not from `final_time`). The time scales linearly with `final_time` and roughly as h^-3.5 with the grid.
- `--long` runs `tests/benchmark_configs_long/*.json` instead (passed on as `BENCH_SET=long`): `zirconium_800` (the production config at 800×800, `final_time` 0.0016), once with each version by default (`--repeat N` for more), time limit 30 min + 30 min per repetition, timeout per run 40 min instead of 10. It shows whether a speedup holds on a large grid with the GPU fully loaded. Its ~10 min on an RTX 2070 SUPER is extrapolated from the 200 and 400 runs (the cost per simulated time grew 9.5× from 100 to 200 and 11.5× from 200 to 400; 13–16× assumed for 800), not yet measured.
- Unlike `test`, `bench` runs when both versions have the same code: that shows the noise. Times depend on the GPU type, so compare speedups, or use `--gpu` to benchmark on the same type.

Plotting (it runs from the repo root and needs `results/<name>/info/setup/config.json`):
```bash
python3 python/basic_plot.py --name MyRun --type both   # selected frames + animations (mp4, or gif without ffmpeg) into results/MyRun/info/plots/
python3 python/basic_plot.py --name MyRun --frames 0 50 last --format gif
python3 python/plot.py --name MyRun --frame 10 --type phase   # --frame last works too, --no-save opens a window instead of saving
```
`python/graphs/*.py` are standalone scripts for thesis figures. They use LaTeX text rendering.

## Architecture

**Data flow:** `config/config.json` → `Parameters::load` → `Model` + `InitialCondition` → TNL `ODESolver` (adaptive Kutta–Merson) loop in `src/main.cu` → `write_frame` (`src/ResultIO.cu`) → `results/<FOLDER>/calculations/Result-NNNNN.npy` → Python plotters.

**State vector:** there is a single TNL `Vector` `u` of length `2*sizeX*sizeY` that lives on the CUDA device. Phase values are in `[0, N)` and concentration values in `[N, 2N)`, with `N = sizeX*sizeY` and linear index `j*sizeX + i`. Use `Model::phase_at` / `conc_at` to read them.

**Grid (`include/Grid.hpp`):** a plain struct with `sizeX`, `sizeY`, `domain` and the steps `hx`, `hy` (computed once by its constructor), plus `points()` (= N), `index(i, j)` and the coordinates `x(i)`, `y(j)`. It is built from the config by `Parameters::grid()`, held by `Model` and `InitialCondition` and passed to `write_frame`/`read_frame` and the IC setters. Use it instead of writing the index arithmetic or computing h again.

**Time stepping (`main.cu`):** the RHS lambda is evaluated with `parallelFor` on interior points only. Boundary conditions are separate `parallelFor`s that set `fu = 0` on the boundary and write Dirichlet or Neumann values **directly into `u`**. The outer loop advances the solver one output frame at a time (`timeStep = (final-initial)/frame_num`) and writes one file per frame. Unless `custom_integration_time_step` is set, the integration step defaults to `min(hx,hy)^2/5`.

**GPU constraints:** `Model` and `InitialCondition` are captured **by value** in device lambdas. Any method called from the RHS, BCs or IC must be `__cuda_callable__`, and their members end up on the device, so they hold only plain data: `Grid`, `ModelParams` (`include/ModelParams.hpp`, the parameters the GPU needs, built by `Parameters::model_params()`) and numbers. `Parameters` itself (strings, run settings) stays on the host. `Model` does not change after construction, so its methods are `const`. Lambdas use `--extended-lambda`, and the build uses relocatable device code (`-rdc=true`).

**Model (`include/Model.hpp`, `src/Model.cu`):**
- Phase RHS: `(1/alpha) * (div T0(grad p) + f_0/ksi^2 - (b/ksi) * p^2(p-1)^2 * F)`, where `T0` contains the anisotropy `psi(theta) = 1 + A cos(m(theta - theta_0))` (`Anisotropy` in `ModelParams.hpp`, shared with the Wulff IC).
- `F` (the force term) is selected by `force_term.type`. `constant` and `radial` are verification cases with analytic solutions (a shrinking circle or Wulff shape). `zirconium` uses the Gibbs-energy difference `G_m_alpha - G_m_beta` from `include/constants.hpp`.
- Concentration only evolves when the force term is `zirconium`; otherwise its RHS is 0.
- Physical constants and thermodynamic functions for Zr–Nb (CALPHAD-style Gibbs energies, diffusivities, mobilities) are in `include/constants.hpp`.

**Config validation (`Parameters::load`):** it requires `m >= 2` and the convexity condition `A <= 1/(m^2-1)`. `radius_proportion` is scaled by the shorter side of the domain.

**Output layout for `results/<FOLDER>/`:**
- `calculations/Result-NNNNN.npy`: the only output format. Array of shape `(2, sizeY, sizeX)` of `Real` (phase, then concentration; x and y are not stored, they are computed from the config). Written and read by `write_frame` / `read_frame` in `src/ResultIO.cu` (using `include/npy.hpp`), read in Python by `python/dataset.py`. (A text format `x y phase conc` existed until 2026-11-11 and was removed.)
- `info/setup/`: written by the binary at the start, through `RunFolder` (`src/RunFolder.cu`), which holds the C++ side of the folder layout (`calculations/`, `info/setup/`, `parameters.txt`) and writes the run information.
  - `parameters.txt`: the parameters, then the run date (start time), computer (host name, e.g. `gp11`) and GPU (written by the binary at the start), the computation time (written at the end) and the git commit (appended by `make run`).
  - `config.json`: a copy of the config, which the Python plotters read. `Dataset` falls back to `info/config.json` for results computed before the subfolders existed.
  - `parameters.tex`
  - `sources.tar.gz`: only for cluster runs, the exact sources the job compiled (with uncommitted changes and the config).
  - `mode.txt`: only for cluster jobs, `run`, `test` or `bench` (written by `run.sbatch`; tests and benchmarks also keep `baseline.tar.gz`).
- `info/plots/`: frames (`.jpg`) and animations (`.mp4`/`.gif`) from the Python plotters.
- `info/boundary/`: `interface_num.crv` and `interface_anal.crv` from `write_boundary_files`.
- The paths of the three subfolders are defined once in `Dataset` (`setup_folder`, `plots_folder`, `boundary_folder`); on the C++ side the paths are in `RunFolder`. The Makefile (`PARAM_PATH`), `tests/common.sh` and `tests/benchmark.sh` also know where `parameters.txt` is.

**Python (`python/plotter.py`):** `Plotter` base class with `BoundaryPlotter2D` (draws the 0.5-level contour, which can be overlaid with the analytic circle or Wulff boundary), `SurfacePlotter` (3D) and `CutPlotter` (1D cut). The boundary computation does not live in the plotters: `python/boundary.py` has the class `Boundary` (numerical curves from contourpy, analytic curves) and `write_boundary_files`, which `basic_plot.py` calls once to write `info/boundary/interface_num.crv` and `interface_anal.crv` for all frames (the files are always rewritten from scratch; `plot.py` never writes them). `python/cli.py` is shared by `basic_plot.py` and `plot.py`: it holds the common arguments, `decide_drawing` (whether to draw concentration plots and the analytic boundary, from `force_term.type` and the `initial_condition` name) and `build_plotters`.

## Common extension points

- **New initial condition:**
  1. Add it to `ICType` in `include/enums.hpp`, before `Count`.
  2. Add both string conversions in `src/enums.cu`.
  3. Add a `case` in `InitialCondition::get_phase` (and `get_concentration` if needed) in `src/InitialCondition.cu`.
  4. If it has an analytic solution, add its name to `ANALYTIC_INITIAL_CONDITIONS` and handle it in `analytic_polyline` (both in `python/boundary.py`; `decide_drawing` in `python/cli.py` uses the list).
- **New config parameter:** add a field in `include/Parameters.hpp` (floating-point fields are `Real`, not `double`, so that a `float` build does no double arithmetic on the GPU), then read it in `Parameters::load`. Also add it to `save_human_readable` and `save_for_latex` so runs stay self-documenting. If the GPU code needs it, add it to `ModelParams` and copy it in `Parameters::model_params()` (or, for an IC, to the members of `InitialCondition`).
- **Initial condition from file:** with `init_cond_from_file: true`, the path is resolved relative to the run's `calculations/` folder and the file must be a `Result-*.npy` frame of float or double (converted to `Real`; the shape must be `(2, sizeY, sizeX)`). `scaling: true` embeds every second point of it (a grid with half the resolution) into the centre of the domain. The code is `set_init_cond_from_file` in `src/InitialCondition.cu`.
