#!/bin/bash
# Submit from the repo root:
#   sbatch --export=ALL,FOLDER=MyRun slurm/run.sbatch
#   sbatch --export=ALL,FOLDER=MyRun,DELETE_CALCS=1 slurm/run.sbatch
#
#SBATCH --job-name=phase-field
#SBATCH --partition=gp               # any GPU node; the Makefile compiles for whichever GPU it lands on
#SBATCH --gres=gpu:1                 # one GPU is enough, the solver is single-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G                     # host RAM; cluster default is 2G (max 60G), too tight for nvcc + matplotlib
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

# Slurm starts the job in the directory you ran sbatch from, which has to be
# the repo root (config/config.json and results/ are relative paths).
cd "${SLURM_SUBMIT_DIR}"

: "${FOLDER:?Set FOLDER, e.g. sbatch --export=ALL,FOLDER=MyRun slurm/run.sbatch}"

# CUDA is installed in /opt/cuda on the cluster (there is no module system).
export CUDA_HOME=/opt/cuda
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
nvcc --version | tail -n 2

nvidia-smi --query-gpu=name,compute_cap --format=csv
# GIT_COMMIT is passed in by cluster.sh (there is no .git folder on the cluster)
extra=()
if [ -n "${GIT_COMMIT:-}" ]; then
    extra+=(GIT_COMMIT="${GIT_COMMIT}")
fi
make run FOLDER="${FOLDER}" CONFIG="${CONFIG:-config/config.json}" DELETE_CALCS="${DELETE_CALCS:-0}" "${extra[@]}"
