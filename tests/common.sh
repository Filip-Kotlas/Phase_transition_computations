# Functions shared by tests/regression.sh and tests/benchmark.sh. Sourced, not run.
# The script that sources it sets $out (its output folder) first.

# Folder of the test scripts (this file's folder), for prepare_input.py, and the source tree it
# belongs to: the harness, which supplies the configs and the comparison scripts. It is the
# version under test, unless that is another commit (cluster.sh test/bench --new COMMIT).
tests_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
harness_dir=$(cd "$tests_dir/.." && pwd)

# A test case that runs longer than this is stuck; the solver loop has hung before (job 53761)
# and would otherwise hold the GPU until the job's time limit.
CASE_TIMEOUT=${CASE_TIMEOUT:-600}

# Progress for `cluster.sh status`: one line in $out/progress, rewritten at every step:
#   build VERSION
#   run K TOTAL START CASE VERSION   (K-th of TOTAL runs, started at START in seconds since 1970)
#   compare
# The calling script deletes it when it ends normally.
progress() { echo "$*" > "$out/progress"; }

# Builds the version in $1 and prints the path of its binary. Only plain `make` is used, which
# every version of the Makefile has (the binary moved from build/ to build/sm_XX/ over time).
build() {
    echo ">> Building $1" >&2
    # Explicit check: set -e does not apply inside $(...), from which this function is called.
    make -C "$1" -j 2 >&2 || { echo "The build of $1 failed" >&2; exit 1; }
    local binary
    binary=$(find "$1/build" -name Phase-field -type f -perm -u+x | head -n 1)
    [ -n "$binary" ] || { echo "No Phase-field binary in $1/build after the build" >&2; exit 1; }
    echo "$binary"
}

# run_case BINARY SRC_DIR CONFIG RUN_DIR COMMIT
# Runs the solver BINARY with CONFIG into the fresh folder RUN_DIR, from SRC_DIR (like `make run`
# does); the output goes to RUN_DIR/solver.log. Sets
#   run_status   exit code of the solver (124 = killed after CASE_TIMEOUT seconds)
#   run_seconds  wall time of the whole solver process in seconds (start-up included)
# and prints the end of solver.log when the run failed. It never fails itself.
run_case() {
    local binary=$1 src=$2 config=$3 run_dir=$4 commit=$5 start end
    rm -rf "$run_dir"
    mkdir -p "$run_dir"
    # The input file of a case that starts from a file (nothing for the other cases).
    python3 "$tests_dir/prepare_input.py" "$config" "$run_dir"

    run_status=0
    start=$(date +%s.%N)
    (cd "$src" && timeout "$CASE_TIMEOUT" "$binary" "$run_dir" "$config") > "$run_dir/solver.log" 2>&1 || run_status=$?
    end=$(date +%s.%N)
    run_seconds=$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.2f", e - s }')

    # parameters.txt gets the commit like with `make run`, so compare_runs.py can show it.
    # (Versions before commit 98e885c write it into info/ instead of info/setup/.)
    local parameters
    for parameters in "$run_dir/info/setup/parameters.txt" "$run_dir/info/parameters.txt"; do
        if [ -f "$parameters" ]; then
            echo "Git commit: $commit" >> "$parameters"
            break
        fi
    done
    if [ "$run_status" != 0 ]; then
        if [ "$run_status" = 124 ]; then echo "   timed out after $CASE_TIMEOUT s"
        else echo "   failed with exit code $run_status"; fi
        echo "   last lines of $run_dir/solver.log:"
        tail -n 5 "$run_dir/solver.log" | sed 's/^/     /'
    fi
    return 0
}
