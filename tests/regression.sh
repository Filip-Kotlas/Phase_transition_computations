#!/bin/bash
# Regression test: builds two versions of the code, runs every tests/configs/*.json with both
# and compares the results with python/compare_runs.py.
#
#   tests/regression.sh NEW_DIR BASE_DIR OUT_DIR [NEW_COMMIT [BASE_COMMIT]]
#
# NEW_DIR and BASE_DIR are unpacked source trees (the version under test and the baseline).
# The configs, this script and the comparison always come from NEW_DIR, so both versions compute
# exactly the same problems. Results go to OUT_DIR/new/<case> and OUT_DIR/base/<case>, the
# comparison to OUT_DIR/report.txt. Exit code 0 when every case is IDENTICAL or CLOSE.
#
# Normally started by slurm/run.sbatch for `slurm/cluster.sh test`, but it works on any machine
# with a GPU, nvcc and TNL.

set -euo pipefail

[ $# -ge 3 ] || { echo "Usage: $0 NEW_DIR BASE_DIR OUT_DIR [NEW_COMMIT [BASE_COMMIT]]" >&2; exit 2; }
new_dir=$(cd "$1" && pwd)
base_dir=$(cd "$2" && pwd)
mkdir -p "$3"
out=$(cd "$3" && pwd)
new_commit="${4:-unknown}"
base_commit="${5:-unknown}"

# A test case that runs longer than this is stuck (each one takes seconds); the solver loop has
# hung before (job 53761) and would otherwise hold the GPU until the job's time limit.
CASE_TIMEOUT=600

configs=("$new_dir"/tests/configs/*.json)
[ -e "${configs[0]}" ] || { echo "No configs in $new_dir/tests/configs" >&2; exit 2; }

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
new_binary=$(build "$new_dir")
base_binary=$(build "$base_dir")

for config in "${configs[@]}"; do
    case_name=$(basename "$config" .json)
    for version in base new; do
        if [ "$version" = new ]; then binary=$new_binary; src=$new_dir; commit=$new_commit
        else                          binary=$base_binary; src=$base_dir; commit=$base_commit; fi
        run_dir="$out/$version/$case_name"
        rm -rf "$run_dir"
        mkdir -p "$run_dir"
        # The input file of a case that starts from a file (nothing for the other cases).
        python3 "$new_dir/tests/prepare_input.py" "$config" "$run_dir"

        echo ">> $case_name: running the $version version ($commit)"
        # Run from the source tree, like `make run` does; the output goes to solver.log.
        status=0
        (cd "$src" && timeout "$CASE_TIMEOUT" "$binary" "$run_dir" "$config") > "$run_dir/solver.log" 2>&1 || status=$?
        # parameters.txt gets the commit like with `make run`, so compare_runs.py can show it.
        [ -f "$run_dir/info/setup/parameters.txt" ] && echo "Git commit: $commit" >> "$run_dir/info/setup/parameters.txt"
        if [ "$status" != 0 ]; then
            [ "$status" = 124 ] && echo "   timed out after $CASE_TIMEOUT s" || echo "   failed with exit code $status"
            echo "   last lines of $run_dir/solver.log:"
            tail -n 5 "$run_dir/solver.log" | sed 's/^/     /'
        fi
    done
done

echo
echo ">> Comparing base ($base_commit) with new ($new_commit)"
# The exit code of the comparison is the exit code of the test.
status=0
python3 "$new_dir/python/compare_runs.py" "$out/base" "$out/new" > "$out/report.txt" || status=$?
cat "$out/report.txt"
exit "$status"
