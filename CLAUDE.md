# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

CUDA/C++17 phase-field solver for a 2D phase transition (anisotropic Allen–Cahn-type equation for the phase `p`, coupled with diffusion of the Nb concentration `c` in a Zr–Nb alloy). It is research code for a thesis at FNSPE CTU Prague. Python scripts post-process and plot the results. Comments and many messages are in Czech.

## Build & run

The code needs `nvcc` and the header-only **TNL** library (Template Numerical Library), which is expected in `~/.local/include`. Builds only work on a Linux machine with an NVIDIA GPU. The Makefile compiles with `-arch=native` (override with `make CUDA_ARCH=...`) and keeps objects in `build/sm_<compute capability>/`, so build on the machine that will run the binary. Jobs for the GPU cluster are submitted from a local machine with `slurm/cluster.sh` (see `slurm/run.sbatch`); from Windows PowerShell use the launcher `.\slurm\cluster.cmd`, which runs the script with Git's bash and Windows OpenSSH (the key lives in the Windows ssh-agent).

```bash
make                                  # build ./build/Phase-field
make run FOLDER=MyRun                 # build, run with config/config.json -> results/MyRun, then plot
make run FOLDER=MyRun DELETE_CALCS=1  # same, but keep only the first and last Result-*.txt afterwards
make run_all                          # run every comp_configurations/*.json (each one is copied over config/config.json)
make clean
```

- The binary is called as `Phase-field [result_folder [config_file]]`; the config defaults to `config/config.json` (`make run CONFIG=...` overrides it). `slurm/cluster.sh submit` uploads a per-run copy to `config/queued/<FOLDER>.json` on the cluster.
- Run everything from the repo root. The binary resolves its default config through a relative path, and the Python scripts look in `./results/<name>`.
- The Makefile has no header dependency tracking. Run `make clean` after changing any `include/*.hpp`, otherwise stale `.o` files cause hard-to-find bugs (this has happened before).
- `make run_all` overwrites the tracked `config/config.json`. `comp_configurations/` is gitignored, apart from `.gitkeep`.
- There are no tests and no linter.

Plotting (it runs from the repo root and needs `results/<name>/info/config.json`):
```bash
python3 python/basic_plot.py --name MyRun --type both   # selected frames + animations (mp4, or gif without ffmpeg) into results/MyRun/info/
python3 python/basic_plot.py --name MyRun --frames 0 50 last --format gif
python3 python/plot.py --name MyRun --frame 10 --type phase   # --frame last works too, --no-save opens a window instead of saving
```
`python/graphs/*.py` are standalone scripts for thesis figures. They use LaTeX text rendering.

## Architecture

**Data flow:** `config/config.json` → `Parameters::load` → `Problem` + `InitialCondition` → TNL `ODESolver` (adaptive Kutta–Merson) loop in `src/main.cu` → `results/<FOLDER>/calculations/Result-NNNNN.txt` → Python plotters.

**State vector:** there is a single TNL `Vector` `u` of length `2*sizeX*sizeY` that lives on the CUDA device. Phase values are in `[0, N)` and concentration values in `[N, 2N)`, with `N = sizeX*sizeY` and linear index `j*sizeX + i`. Use `Problem::phase_at` / `conc_at` to read them.

**Time stepping (`main.cu`):** the RHS lambda is evaluated with `parallelFor` on interior points only. Boundary conditions are separate `parallelFor`s that set `fu = 0` on the boundary and write Dirichlet or Neumann values **directly into `u`**. The outer loop advances the solver one output frame at a time (`timeStep = (final-initial)/frame_num`) and writes one file per frame. Unless `custom_integration_time_step` is set, the integration step defaults to `min(hx,hy)^2/5`.

**GPU constraints:** `Problem` and `InitialCondition` are captured **by value** in device lambdas. Any method called from the RHS, BCs or IC must be `__cuda_callable__`, and their members (including the `Parameters` copy) end up on the device. Lambdas use `--extended-lambda`, and the build uses relocatable device code (`-rdc=true`).

**Model (`src/Problem.cu`):**
- Phase RHS: `(1/alpha) * (div T0(grad p) + f_0/ksi^2 - (b/ksi) * p^2(p-1)^2 * F)`, where `T0` contains the anisotropy `psi(theta) = 1 + A cos(m(theta - theta_0))`.
- `F` (the force term) is selected by `force_term.type`. `constant` and `radial` are verification cases with analytic solutions (a shrinking circle or Wulff shape). `zirconium` uses the Gibbs-energy difference `G_m_alpha - G_m_beta` from `include/constants.hpp`.
- Concentration only evolves when the force term is `zirconium`; otherwise its RHS is 0.
- Physical constants and thermodynamic functions for Zr–Nb (CALPHAD-style Gibbs energies, diffusivities, mobilities) are in `include/constants.hpp`.
- The `#define`s at the top of `Problem.cu` (`INIT`, `FORCE`, `C_TEST_X`, …) are legacy. Runtime behaviour is controlled by the JSON config.

**Config validation (`Parameters::load`):** it requires `m >= 2` and the convexity condition `A <= 1/(m^2-1)`. `radius_proportion` is scaled by the shorter side of the domain.

**Output layout for `results/<FOLDER>/`:**
- `calculations/Result-NNNNN.txt`: lines of `x y phase conc`, with a blank line after each grid row.
- `info/`:
  - `parameters.txt`: `make run` appends the run date and git commit to it.
  - `config.json`: a copy of the config, which the Python plotters read.
  - `parameters.tex`
  - generated PNGs and GIFs.

**Python (`python/plotter.py`):** `Plotter` base class with `BoundaryPlotter2D` (the 0.5-level contour, which can be overlaid with the analytic circle or Wulff boundary and saved as boundary points), `SurfacePlotter` (3D) and `CutPlotter` (1D cut). `python/cli.py` is shared by `basic_plot.py` and `plot.py`: it holds the common arguments, `decide_drawing` (whether to draw concentration plots and the analytic boundary, from `force_term.type` and the `initial_condition` name) and `build_plotters`.

## Common extension points

- **New initial condition:**
  1. Add it to `ICType` in `include/enums.hpp`, before `Count`.
  2. Add both string conversions in `src/enums.cu`.
  3. Add a `case` in `InitialCondition::get_phase` (and `get_concentration` if needed) in `src/InitialCondition.cu`.
  4. If it has an analytic solution, update `decide_drawing` in `python/cli.py` and `BoundaryPlotter2D`.
- **New config parameter:** add a field in `include/Parameters.hpp`, then read it in `Parameters::load`. Also add it to `save_human_readable` and `save_for_latex` so runs stay self-documenting.
- **Initial condition from file:** with `init_cond_from_file: true`, the path is resolved relative to the run's `calculations/` folder and the file must use the `Result-*.txt` format. `scaling: true` embeds a grid with half the resolution into the centre of the domain.
