#!/bin/bash
# Benchmark: how much faster (or slower) a version of the code is than a baseline. Builds both
# versions, runs every tests/benchmark_configs/*.json REPEAT times with each of them on the same
# GPU and compares the computation times with python/compare_times.py.
#
#   tests/benchmark.sh NEW_DIR BASE_DIR OUT_DIR [NEW_COMMIT [BASE_COMMIT]]
#   REPEAT=5 tests/benchmark.sh ...          (default REPEAT=3)
#   BENCH_SET=long tests/benchmark.sh ...     the long benchmark: tests/benchmark_configs_long/*.json
#                                             (cluster.sh bench --long, which passes REPEAT=1)
#
# NEW_DIR and BASE_DIR are unpacked source trees, as for tests/regression.sh; the configs and the
# comparison come from the tree of this script (normally NEW_DIR itself), each version runs the
# configs of its config set (BASE_CONFIG_SET, NEW_CONFIG_SET), and a case that the base version's
# set does not have is skipped, as in tests/regression.sh. In OUT_DIR:
#   times.csv          one line per run: case, version, repetition, computation time (the solver's
#                      own "Computation time" from parameters.txt: the time loop including frame
#                      writing, without start-up, config and initial condition; empty in versions
#                      that do not write it), wall time of the whole process, exit code
#   base/<case>, new/<case>   the first repetition of each case, complete (frames, parameters.txt,
#                      solver.log). python/compare_runs.py compares them into comparison.txt, which
#                      shows whether the faster version still computes the same thing.
#   repeats/<version>/<case>-<k>   the other repetitions, without their frames
#   report.txt         the timing table
# Exit code 0 when every case has times for both versions.
#
# Noise: the versions alternate (base new, new base, base new, ...), so a slow drift on the node
# (another job, the GPU warming up) hits both alike, and the report uses the median of the
# repetitions. The cases are short versions of real runs: the time per frame is constant during
# a run, so a few frames measure the same thing as a whole run.
#
# Normally started by slurm/run.sbatch for `slurm/cluster.sh bench`.

set -euo pipefail

[ $# -ge 3 ] || { echo "Usage: [REPEAT=N] [BENCH_SET=long] $0 NEW_DIR BASE_DIR OUT_DIR [NEW_COMMIT [BASE_COMMIT]]" >&2; exit 2; }
new_dir=$(cd "$1" && pwd)
base_dir=$(cd "$2" && pwd)
mkdir -p "$3"
out=$(cd "$3" && pwd)
new_commit="${4:-unknown}"
base_commit="${5:-unknown}"
repeat="${REPEAT:-3}"
case "$repeat" in ""|*[!0-9]*|0) echo "REPEAT has to be a positive number, not '$repeat'" >&2; exit 2 ;; esac

# The cases of the long benchmark take about 10 minutes each, longer than the default timeout
# for a stuck case in common.sh.
case "${BENCH_SET:-}" in
    "")   configs_folder=benchmark_configs ;;
    long) configs_folder=benchmark_configs_long; CASE_TIMEOUT=${CASE_TIMEOUT:-2400} ;;
    *)    echo "BENCH_SET has to be empty or 'long', not '$BENCH_SET'" >&2; exit 2 ;;
esac

# progress, build, run_case, config_folder and CASE_TIMEOUT.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

new_configs=$(config_folder "${NEW_CONFIG_SET:-}" "$configs_folder") || exit 2
base_configs=$(config_folder "${BASE_CONFIG_SET:-}" "$configs_folder") || exit 2
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

times="$out/times.csv"
echo "case,version,repetition,computation_s,wall_s,exit_code" > "$times"

total=$(( 2 * repeat * ${#cases[@]} ))
k=0
for case_name in "${cases[@]}"; do
    for (( r = 1; r <= repeat; r++ )); do
        # Odd repetitions run base first, even ones new first.
        if (( r % 2 == 1 )); then order="base new"; else order="new base"; fi
        for version in $order; do
            if [ "$version" = new ]; then binary=$new_binary; src=$new_dir; commit=$new_commit; config=$new_configs/$case_name.json
            else                          binary=$base_binary; src=$base_dir; commit=$base_commit; config=$base_configs/$case_name.json; fi
            k=$(( k + 1 ))
            progress run "$k" "$total" "$(date +%s)" "$case_name" "$version"
            echo ">> [$k/$total] $case_name, repetition $r: running the $version version ($commit)"

            if (( r == 1 )); then run_dir="$out/$version/$case_name"
            else                  run_dir="$out/repeats/$version/$case_name-$r"; fi
            run_case "$binary" "$src" "$config" "$run_dir" "$commit"

            # "Computation time: 0h 00m 29s (29.3 s)" -> 29.3; empty when the run did not write it
            # (|| true: a run that crashed early has no parameters.txt, and set -e would stop here).
            computation=$(sed -n 's/^Computation time:.*(\([0-9.]*\) s)$/\1/p' "$run_dir/info/setup/parameters.txt" 2>/dev/null | tail -n 1 || true)
            echo "$case_name,$version,$r,$computation,$run_seconds,$run_status" >> "$times"
            if [ -n "$computation" ]; then echo "   computation time $computation s, wall time $run_seconds s"
            else                           echo "   computation time not recorded, wall time $run_seconds s"; fi
            # The frames of the first repetition are enough for the comparison of the results.
            (( r == 1 )) || rm -rf "$run_dir/calculations"
        done
    done
done

echo
echo ">> Comparing the results of base ($base_commit) and new ($new_commit)"
progress compare
python3 "$harness_dir/python/compare_runs.py" "$out/base" "$out/new" > "$out/comparison.txt" || true
tail -n 1 "$out/comparison.txt"

echo
echo ">> Comparing the times"
status=0
{
    echo "Configs: base ${base_configs#$harness_dir/}, new ${new_configs#$harness_dir/}"
    [ ${#skipped[@]} -eq 0 ] || echo "Not run, the base version's configs do not have them: ${skipped[*]}"
    echo
    python3 "$harness_dir/python/compare_times.py" "$out" || status=$?
} > "$out/report.txt"
cat "$out/report.txt"
rm -f "$out/progress"
exit "$status"
