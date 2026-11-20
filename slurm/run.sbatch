#!/bin/bash
# Submitted by `slurm/cluster.sh submit`, not by hand. cluster.sh uploads the sources of the run
# as jobs/<FOLDER>/sources.tar.gz into the cluster folder, runs sbatch in that folder and passes
# FOLDER, DELETE_CALCS and GIT_COMMIT (the job name is set to FOLDER).
#
# `slurm/cluster.sh test` submits it with MODE=test: it also uploads the baseline version as
# jobs/<FOLDER>/baseline.tar.gz and passes BASELINE_COMMIT. The job then runs the regression
# test (tests/regression.sh) instead of `make run`. `slurm/cluster.sh bench` is the same with
# MODE=bench (and REPEAT, and BENCH_SET=long for --long): the job runs the benchmark
# (tests/benchmark.sh).
# With --new COMMIT, test and bench also upload that commit as jobs/<FOLDER>/new.tar.gz and pass
# NEW_COMMIT: it is then the version under test, and sources.tar.gz (your working tree) only
# supplies this script, the test scripts, the configs and the comparison.
# The mode is kept in results/<FOLDER>/info/setup/mode.txt, for `cluster.sh status` and `list`.
#
#SBATCH --partition=gp               # any GPU node; the Makefile compiles for whichever GPU it lands on
#SBATCH --gres=gpu:1                 # one GPU is enough, the solver is single-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
# Host RAM: the cluster default 2G (max 60G) is too tight for nvcc + matplotlib. Jobs peaked at 1.7 GB
# (sacct MaxRSS, 2026-11-12); raise it (or use cluster.sh --mem) if a job is killed for lack of memory.
#SBATCH --mem=4G
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
mode="${MODE:-run}"
case "$mode" in
    run)  archives="sources.tar.gz" ;;
    test) archives="sources.tar.gz baseline.tar.gz"
          echo "Regression test: ${GIT_COMMIT:-unknown} against the baseline ${BASELINE_COMMIT:-unknown}" ;;
    bench) archives="sources.tar.gz baseline.tar.gz"
          echo "Benchmark${BENCH_SET:+ ($BENCH_SET)}: ${GIT_COMMIT:-unknown} against the baseline ${BASELINE_COMMIT:-unknown}, ${REPEAT:-3} run(s) of each" ;;
    *)    echo "Unknown MODE '$mode'" >&2; exit 1 ;;
esac
if [ "$mode" != run ] && [ -n "${NEW_COMMIT:-}" ]; then
    archives="$archives new.tar.gz"
    echo "The version under test is the commit $NEW_COMMIT, not the uploaded working tree."
fi
job_dir="$base/jobs/$FOLDER"
uploaded() {
    local a
    for a in $archives; do [ -f "$job_dir/$a" ] || return 1; done
}
# The home folder reaches the nodes over NFS, which caches folder contents for up to about a
# minute. A job that starts right after submit can therefore miss the fresh upload, so wait for
# it (listing the folder makes NFS read it again).
for _ in $(seq 24); do
    ls "$base/jobs" "$job_dir" >/dev/null 2>&1 || true
    uploaded && break
    sleep 5
done
uploaded || { echo "Missing $archives in $job_dir after waiting 2 minutes (was the job submitted by cluster.sh?)" >&2; exit 1; }

# Every job compiles and runs its own copy of the sources, unpacked on the node's local disk.
# So a later submit cannot change the code of a queued job, and jobs cannot disturb each other.
# A test also unpacks the baseline version, into its own folder, and the version under test, if
# it is a commit (otherwise the version under test is the working tree in $work).
work="${TMPDIR:-/tmp}/phase-field-${SLURM_JOB_ID}"
base_work="$work-baseline"
new_work="$work"
trap 'rm -rf "$work" "$base_work" "$work-new"' EXIT
trap 'exit 143' TERM                 # scancel or time limit: exit, so that the EXIT cleanup runs
mkdir -p "$work"
tar xzf "$job_dir/sources.tar.gz" -C "$work"
if [ "$mode" != run ]; then
    mkdir -p "$base_work"
    tar xzf "$job_dir/baseline.tar.gz" -C "$base_work"
    if [ -n "${NEW_COMMIT:-}" ]; then
        new_work="$work-new"
        mkdir -p "$new_work"
        tar xzf "$job_dir/new.tar.gz" -C "$new_work"
    fi
fi

# Keep the exact sources with the results. This also creates results/$FOLDER right away, so the
# run name stays taken. The upload folder is not needed any more after that.
mkdir -p "results/$FOLDER/info/setup"
for a in $archives; do
    cp "$job_dir/$a" "results/$FOLDER/info/setup/$a"
done
echo "$mode" > "results/$FOLDER/info/setup/mode.txt"
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
# Every mode takes its binaries from the build cache in the cluster folder when the same code was
# already built for this GPU type, and stores the ones it builds there (slurm/cached_build.sh;
# test and bench through build() in tests/common.sh). `cluster.sh cache` lists and prunes it.
export BUILD_CACHE="$base/build_cache"
if [ "$mode" = test ]; then
    # Both versions run on this GPU, so differences between GPU types cannot show up as
    # differences between the versions. The results and report.txt go to results/$FOLDER.
    bash tests/regression.sh "$new_work" "$base_work" "results/$FOLDER" "${GIT_COMMIT:-unknown}" "${BASELINE_COMMIT:-unknown}"
elif [ "$mode" = bench ]; then
    # The same for the times: both versions on this GPU, report.txt in results/$FOLDER.
    REPEAT="${REPEAT:-3}" BENCH_SET="${BENCH_SET:-}" bash tests/benchmark.sh "$new_work" "$base_work" "results/$FOLDER" "${GIT_COMMIT:-unknown}" "${BASELINE_COMMIT:-unknown}"
else
    # E.g. the other runs of a submit-many take the binary from the cache.
    binary=$(bash slurm/cached_build.sh "$work" "$BUILD_CACHE" "${GIT_COMMIT:-unknown}")
    # The copy has no .git folder, so the commit comes from cluster.sh. config/config.json is the
    # config from the uploaded sources, i.e. the one that was in your repo when you submitted
    # (or the file given with submit --config / submit-many, which cluster.sh put there).
    make run FOLDER="${FOLDER}" BINARY="$binary" CONFIG=config/config.json DELETE_CALCS="${DELETE_CALCS:-0}" GIT_COMMIT="${GIT_COMMIT:-unknown}"
fi
