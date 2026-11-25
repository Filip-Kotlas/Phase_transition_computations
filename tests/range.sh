#!/bin/bash
# Range test: which commit of a range changed the results? Builds every version of a chain of
# commits C0, C1, ..., Cn (each the parent of the next), runs the regression test cases with each
# and compares every commit with its parent (python/compare_runs.py).
#
#   tests/range.sh OUT_DIR DIR0 COMMIT0 SET0 DIR1 COMMIT1 SET1 [DIR COMMIT SET]...
#
# DIRi is the unpacked source tree of COMMITi (a short hash, also the folder name of its runs) and
# SETi its config set in tests/config_sets/ ("-": the current tests/configs; cluster.sh chooses it
# with config_set_of). C0 is only the baseline of C1. The configs and the comparison come from the
# tree of this script (the harness: the uploaded working tree). With RANGE_SUBJECTS set to a file
# of lines "COMMIT SUBJECT", the report shows the subjects of the commits.
#
# Each version runs once, every case of its config set, into OUT_DIR/runs/COMMIT/<case>: its runs
# are the new version of its own pair and the baseline of the next one. A commit whose binary is
# bitwise the same as its parent's, with the same config set, is not run (its results would be the
# same: the solver is deterministic on one GPU); its verdict is SAME BINARY, and the next commit is
# compared with the parent's runs. A version that does not build is not run, and its pairs are NOT
# COMPARED. A pair compares the cases that both config sets have (the report names the others,
# e.g. hydride, which the sets of the commits before it do not have). Each comparison is written to
# OUT_DIR/pairs/COMMIT.txt as soon as both versions have run (so a job that hits its time limit
# keeps them), everything to OUT_DIR/report.txt at the end: one entry per commit, then the
# comparisons. Its last line names the commits that changed the results (cluster.sh list shows it).
# Exit code 0 when every commit is IDENTICAL, CLOSE or SAME BINARY against its parent.
#
# Normally started by slurm/run.sbatch for `slurm/cluster.sh test NAME --range A..B`, but it works
# on any machine with a GPU, nvcc and TNL.

set -euo pipefail

usage() { echo "Usage: $0 OUT_DIR DIR0 COMMIT0 SET0 DIR1 COMMIT1 SET1 [DIR COMMIT SET]..." >&2; exit 2; }
[ $# -ge 7 ] && [ $(( ($# - 1) % 3 )) -eq 0 ] || usage
mkdir -p "$1"
out=$(cd "$1" && pwd)
shift

# progress, build, run_case, config_folder and CASE_TIMEOUT.
source "$(dirname "${BASH_SOURCE[0]}")/common.sh"

dirs=(); commits=(); sets=(); configs=()
while [ $# -gt 0 ]; do
    dir=$(cd "$1" && pwd) || exit 2
    set_name=$3; [ "$set_name" = - ] && set_name=""
    folder=$(config_folder "$set_name" configs) || exit 2
    dirs+=("$dir"); commits+=("$2"); sets+=("$set_name"); configs+=("$folder")
    shift 3
done
n=${#commits[@]}

subject_of() {
    [ -n "${RANGE_SUBJECTS:-}" ] && [ -f "$RANGE_SUBJECTS" ] || return 0
    sed -n "s/^$1 //p" "$RANGE_SUBJECTS" | head -n 1
}

# The case names of a config folder.
cases_in() {
    local config
    for config in "$1"/*.json; do
        [ -e "$config" ] && basename "$config" .json
    done
    return 0
}

# Upper bound of the number of runs for the progress line; lowered when a version is not run.
total=0
for folder in "${configs[@]}"; do
    total=$(( total + $(cases_in "$folder" | wc -l) ))
done
k=0

binaries=(); runs_dirs=(); verdicts=(); details=()
mkdir -p "$out/pairs"

# compare_with_parent I: compares version I with version I-1 and sets verdicts[I] and details[I].
#   verdicts: unchanged, changed, same (SAME BINARY), not_compared
compare_with_parent() {
    local i=$1 p=$(( $1 - 1 )) case_name status summary
    local commit=${commits[$i]} parent=${commits[$p]}
    local pair_file="$out/pairs/$commit.txt"
    # The cases both config sets have; the others are named in the report.
    local both=() skipped=()
    for case_name in $(cases_in "${configs[$i]}"); do
        if [ -f "${configs[$p]}/$case_name.json" ]; then both+=("$case_name"); else skipped+=("$case_name"); fi
    done
    for case_name in $(cases_in "${configs[$p]}"); do
        [ -f "${configs[$i]}/$case_name.json" ] || skipped+=("$case_name")
    done
    if [ ${#both[@]} -eq 0 ]; then
        verdicts[$i]=not_compared
        details[$i]="NOT COMPARED: the config sets of $commit and $parent have no case in common"
        echo "${details[$i]}" > "$pair_file"
        return 0
    fi
    progress compare
    echo ">> Comparing $commit with its parent $parent"
    status=0
    {
        echo "$commit against its parent $parent"
        echo "Configs: base ${configs[$p]#$harness_dir/}, new ${configs[$i]#$harness_dir/}"
        [ ${#skipped[@]} -eq 0 ] || echo "Not compared, only one of the two config sets has them: ${skipped[*]}"
        echo
        python3 "$harness_dir/python/compare_runs.py" "${runs_dirs[$p]}" "${runs_dirs[$i]}" --cases "${both[@]}" 2>&1 \
            || status=$?
    } > "$pair_file"
    # "4 IDENTICAL, 2 DIFFERENT (zirconium, from_file)", from the "=== case ===" and "Verdict: X"
    # lines that compare_runs.py prints for every case.
    summary=$(awk '
        /^=== .* ===$/ { name = substr($0, 5, length($0) - 8) }
        /^Verdict: /   { v = $2; count[v]++; names[v] = names[v] (names[v] == "" ? "" : ", ") name }
        END {
            split("IDENTICAL CLOSE DIFFERENT FAILED INCOMPATIBLE", order, " ")
            for (j = 1; j <= 5; j++) {
                v = order[j]
                if (!(v in count)) continue
                s = s (s == "" ? "" : ", ") count[v] " " v (v == "IDENTICAL" ? "" : " (" names[v] ")")
            }
            print s
        }' "$pair_file")
    if [ -z "$summary" ] || [ "$status" -gt 1 ]; then
        verdicts[$i]=not_compared
        details[$i]="NOT COMPARED: compare_runs.py failed (exit code $status), see pairs/$commit.txt"
    elif [ "$status" = 0 ]; then
        verdicts[$i]=unchanged; details[$i]=$summary
    else
        verdicts[$i]=changed; details[$i]=$summary
    fi
    [ ${#skipped[@]} -eq 0 ] || details[$i]="${details[$i]}; not compared (in one config set only): ${skipped[*]}"
    return 0
}

for i in "${!commits[@]}"; do
    commit=${commits[$i]}
    p=$(( i - 1 ))
    echo
    echo ">> Version $(( i + 1 ))/$n: $commit $(subject_of "$commit")"
    progress build "$commit"
    if binary=$(build "${dirs[$i]}" "$commit"); then
        binaries[$i]=$binary
    else
        binaries[$i]=""
        echo "   the build of $commit failed, so it is not run"
    fi

    runs_dirs[$i]="$out/runs/$commit"
    same=0
    if [ "$i" -gt 0 ] && [ -n "${binaries[$i]}" ] && [ -n "${binaries[$p]}" ] \
        && [ "${sets[$i]}" = "${sets[$p]}" ] && cmp -s "${binaries[$i]}" "${binaries[$p]}"; then
        same=1
        # The next commit is compared with the runs of the parent.
        runs_dirs[$i]=${runs_dirs[$p]}
    fi

    if [ -z "${binaries[$i]}" ] || [ "$same" = 1 ]; then
        total=$(( total - $(cases_in "${configs[$i]}" | wc -l) ))
        [ "$same" = 0 ] || echo "   the same binary and configs as its parent ${commits[$p]}, so it is not run again"
    else
        for case_name in $(cases_in "${configs[$i]}"); do
            k=$(( k + 1 ))
            progress run "$k" "$total" "$(date +%s)" "$case_name" "$commit"
            echo ">> [$k/$total] $case_name: running $commit"
            run_case "${binaries[$i]}" "${dirs[$i]}" "${configs[$i]}/$case_name.json" "${runs_dirs[$i]}/$case_name" "$commit"
        done
    fi

    [ "$i" -gt 0 ] || continue
    if [ -z "${binaries[$i]}" ] || [ -z "${binaries[$p]}" ]; then
        verdicts[$i]=not_compared
        failed=$commit; [ -n "${binaries[$p]}" ] || failed=${commits[$p]}
        details[$i]="NOT COMPARED: the build of $failed failed"
        echo "${details[$i]}" > "$out/pairs/$commit.txt"
    elif [ "$same" = 1 ]; then
        verdicts[$i]=same
        details[$i]="SAME BINARY as ${commits[$p]} (and the same configs), not run again"
        echo "${details[$i]}" > "$out/pairs/$commit.txt"
    else
        compare_with_parent "$i"
    fi
done

# The report: one entry per commit (subject, then the result against its parent), then the full
# comparisons, and a last line with the commits that changed the results.
changed=(); not_compared=(); same_count=0
for i in $(seq 1 $(( n - 1 ))); do
    case "${verdicts[$i]}" in
        changed)      changed+=("${commits[$i]}") ;;
        not_compared) not_compared+=("${commits[$i]}") ;;
        same)         same_count=$(( same_count + 1 )) ;;
    esac
done
{
    echo "Range test: $(( n - 1 )) commit(s), each against its parent (${commits[0]}..${commits[$(( n - 1 ))]})"
    echo
    for i in "${!commits[@]}"; do
        echo "${commits[$i]}  $(subject_of "${commits[$i]}")"
        if [ "$i" = 0 ]; then line="baseline of the range"; else line=${details[$i]}; fi
        [ -z "${sets[$i]}" ] || line="$line   [configs: tests/config_sets/${sets[$i]}]"
        echo "         $line"
    done
    for i in $(seq 1 $(( n - 1 ))); do
        echo
        echo "===================== ${commits[$i]} against ${commits[$(( i - 1 ))]} ====================="
        cat "$out/pairs/${commits[$i]}.txt"
    done
    echo
    result="no commit changed the results"
    [ ${#changed[@]} -eq 0 ] || result="${#changed[@]} of $(( n - 1 )) commit(s) changed the results: ${changed[*]}"
    [ ${#not_compared[@]} -eq 0 ] || result="$result; not compared: ${not_compared[*]}"
    [ "$same_count" = 0 ] || result="$result ($same_count with the same binary as the parent)"
    echo "Overall: $result"
} > "$out/report.txt"
echo
cat "$out/report.txt"
rm -f "$out/progress"
[ ${#changed[@]} -eq 0 ] && [ ${#not_compared[@]} -eq 0 ]
