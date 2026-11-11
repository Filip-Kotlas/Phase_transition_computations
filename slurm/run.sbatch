#!/bin/bash
# Submitted by `slurm/cluster.sh submit`, not by hand. cluster.sh uploads the sources of the run
# as jobs/<FOLDER>/sources.tar.gz into the cluster folder, runs sbatch in that folder and passes
# FOLDER, DELETE_CALCS and GIT_COMMIT (the job name is set to FOLDER).
#
#SBATCH --partition=gp               # any GPU node; the Makefile compiles for whichever GPU it lands on
#SBATCH --gres=gpu:1                 # one GPU is enough, the solver is single-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G                     # host RAM; cluster default is 2G (max 60G), too tight for nvcc + matplotlib
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%j.out

set -euo pipefail

# Slurm starts the job in the directory sbatch was run from: the cluster folder, which holds
# results/ and jobs/.
cd "${SLURM_SUBMIT_DIR}"
base="$(pwd)"

: "${FOLDER:?FOLDER is not set; submit with slurm/cluster.sh submit}"
# `cluster.sh list` finds the log of a run by this line.
echo "Run $FOLDER, job ${SLURM_JOB_ID} on ${SLURMD_NODENAME:-unknown node}"
job_dir="$base/jobs/$FOLDER"
# The home folder reaches the nodes over NFS, which caches folder contents for up to about a
# minute. A job that starts right after submit can therefore miss the fresh upload, so wait for
# it (listing the folder makes NFS read it again).
for _ in $(seq 24); do
    ls "$base/jobs" "$job_dir" >/dev/null 2>&1 || true
    [ -f "$job_dir/sources.tar.gz" ] && break
    sleep 5
done
[ -f "$job_dir/sources.tar.gz" ] || { echo "Missing $job_dir/sources.tar.gz after waiting 2 minutes (was the job submitted by cluster.sh?)" >&2; exit 1; }

# Every job compiles and runs its own copy of the sources, unpacked on the node's local disk.
# So a later submit cannot change the code of a queued job, and jobs cannot disturb each other.
work="${TMPDIR:-/tmp}/phase-field-${SLURM_JOB_ID}"
trap 'rm -rf "$work"' EXIT
trap 'exit 143' TERM                 # scancel or time limit: exit, so that the EXIT cleanup runs
mkdir -p "$work"
tar xzf "$job_dir/sources.tar.gz" -C "$work"

# Keep the exact sources with the results. This also creates results/$FOLDER right away, so the
# run name stays taken. The upload folder is not needed any more after that.
mkdir -p "results/$FOLDER/info/setup"
cp "$job_dir/sources.tar.gz" "results/$FOLDER/info/setup/sources.tar.gz"
rm -rf "$job_dir"

# The binary and the Python scripts write to ./results, so point it to the shared results folder.
ln -s "$base/results" "$work/results"
cd "$work"

# CUDA is installed in /opt/cuda on the cluster (there is no module system).
export CUDA_HOME=/opt/cuda
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
nvcc --version | tail -n 2

nvidia-smi --query-gpu=name,compute_cap --format=csv
# The copy has no .git folder, so the commit comes from cluster.sh. config/config.json is the
# config from the uploaded sources, i.e. the one that was in your repo when you submitted.
make run FOLDER="${FOLDER}" CONFIG=config/config.json DELETE_CALCS="${DELETE_CALCS:-0}" GIT_COMMIT="${GIT_COMMIT:-unknown}"
