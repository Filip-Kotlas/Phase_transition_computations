#!/bin/bash
# Regression test: builds two versions of the code, runs every tests/configs/*.json with both
# and compares the results with python/compare_runs.py.
#
#   tests/regression.sh NEW_DIR BASE_DIR OUT_DIR [NEW_COMMIT [BASE_COMMIT]]
#
# NEW_DIR and BASE_DIR are unpacked source trees (the version under test and the baseline).
# The configs and the comparison always come from the tree of this script (normally NEW_DIR
# itself; another tree when the version under test is an older commit without them). Each version
# runs the configs of its config set (BASE_CONFIG_SET, NEW_CONFIG_SET; empty: the current
# tests/configs), which hold the same problems in the layout that version reads. A case that the
# base version's set does not have is skipped. Results go to OUT_DIR/new/<case> and
# OUT_DIR/base/<case>, the comparison to OUT_DIR/report.txt. Exit code 0 when every case is
# IDENTICAL or CLOSE.
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

# progress, build, run_case, config_folder and CASE_TIMEOUT (each case takes seconds, so 600 s
# means it is stuck).
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

new_configs=$(config_folder "${NEW_CONFIG_SET:-}" configs) || exit 2
base_configs=$(config_folder "${BASE_CONFIG_SET:-}" configs) || exit 2
# The cases: the configs of the new version that the base version's set has too.
cases=(); skipped=()
for config in "$new_configs"/*.json; do
    [ -e "$config" ] || continue
    case_name=$(basename "$config" .json)
    if [ -f "$base_configs/$case_name.json" ]; then cases+=("$case_name"); else skipped+=("$case_name"); fi
done
[ ${#cases[@]} -gt 0 ] || { echo "No configs in both $new_configs and $base_configs" >&2; exit 2; }
echo ">> Configs: base ${base_configs#$harness_dir/}, new ${new_configs#$harness_dir/}"
[ ${#skipped[@]} -eq 0 ] || echo ">> Skipped, the base version's configs do not have them: ${skipped[*]}"

progress build new
new_binary=$(build "$new_dir" "$new_commit")
progress build base
base_binary=$(build "$base_dir" "$base_commit")

total=$(( 2 * ${#cases[@]} ))
k=0

for case_name in "${cases[@]}"; do
    for version in base new; do
        if [ "$version" = new ]; then binary=$new_binary; src=$new_dir; commit=$new_commit; config=$new_configs/$case_name.json
        else                          binary=$base_binary; src=$base_dir; commit=$base_commit; config=$base_configs/$case_name.json; fi
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
{
    echo "Configs: base ${base_configs#$harness_dir/}, new ${new_configs#$harness_dir/}"
    [ ${#skipped[@]} -eq 0 ] || echo "Not run, the base version's configs do not have them: ${skipped[*]}"
    echo
    python3 "$harness_dir/python/compare_runs.py" "$out/base" "$out/new" || status=$?
} > "$out/report.txt"
cat "$out/report.txt"
rm -f "$out/progress"
exit "$status"
