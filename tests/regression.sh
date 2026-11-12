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

# progress, build, run_case and CASE_TIMEOUT (each case takes seconds, so 600 s means it is stuck).
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

configs=("$new_dir"/tests/configs/*.json)
[ -e "${configs[0]}" ] || { echo "No configs in $new_dir/tests/configs" >&2; exit 2; }

progress build new
new_binary=$(build "$new_dir")
progress build base
base_binary=$(build "$base_dir")

total=$(( 2 * ${#configs[@]} ))
k=0

for config in "${configs[@]}"; do
    case_name=$(basename "$config" .json)
    for version in base new; do
        if [ "$version" = new ]; then binary=$new_binary; src=$new_dir; commit=$new_commit
        else                          binary=$base_binary; src=$base_dir; commit=$base_commit; fi
        k=$(( k + 1 ))
        progress run "$k" "$total" "$(date +%s)" "$case_name" "$version"
        echo ">> [$k/$total] $case_name: running the $version version ($commit)"
        run_case "$binary" "$src" "$config" "$out/$version/$case_name" "$commit"
    done
done

echo
echo ">> Comparing base ($base_commit) with new ($new_commit)"
# The exit code of the comparison is the exit code of the test.
progress compare
status=0
python3 "$new_dir/python/compare_runs.py" "$out/base" "$out/new" > "$out/report.txt" || status=$?
cat "$out/report.txt"
rm -f "$out/progress"
exit "$status"
