#!/bin/bash
# Run from your own computer (Git Bash / WSL / Linux), from anywhere in the repo.
#
#   slurm/cluster.sh submit MyRun [--delete-calcs] [--force] [--gpu TYPE]
#                                                              upload code, sbatch the job
#                                                              (--gpu rtx5080: only that card type)
#   slurm/cluster.sh test MyTest [--against COMMIT] [--force] [--gpu TYPE]
#                                                              regression test: your working tree
#                                                              against COMMIT (default HEAD)
#   slurm/cluster.sh status                                    show your queue
#   slurm/cluster.sh list                                      all runs on the cluster and their state
#   slurm/cluster.sh log                                       tail the newest slurm-*.out
#   slurm/cluster.sh gpus                                      GPU types and how many are free now
#   slurm/cluster.sh fetch [MyRun] [--no-calcs] [--force]      download results/MyRun, or all runs
#   slurm/cluster.sh cancel JOBID
#   slurm/cluster.sh clean                                     delete leftover upload folders jobs/*
#   slurm/cluster.sh delete-results [--yes]                    delete ALL results on the cluster
#
# Every run gets its own copy of the code: submit uploads the sources (including
# config/config.json) as jobs/MyRun/sources.tar.gz, and the job compiles and runs that copy
# on the node's local disk. A later submit therefore never changes the code of a queued job.
# The job keeps the archive as results/MyRun/info/setup/sources.tar.gz and deletes jobs/MyRun.
# submit refuses a name that is already used on the cluster (finished or still queued);
# --force deletes the old results and upload of that name first. It does NOT stop a job
# that is still running under that name, so check `status` first.
# test uploads the working tree like submit, plus the committed version COMMIT (git archive) as
# the baseline. The job builds both, runs every tests/configs/*.json with both on the same GPU
# and compares them (tests/regression.sh, python/compare_runs.py); `log` shows the result and
# results/MyTest/report.txt keeps it. So test your changes BEFORE committing them: once they are
# committed, HEAD is the new code itself, and test refuses to compare the code with itself
# (use --against HEAD~1, or the commit before your changes, then).
# list shows each run in results/ and jobs/ with its state: queued, running (with progress),
# finished (computation time; warns when the animations are missing) or stopped (no job any
# more and not finished: crashed, cancelled or time limit; it names the log to look at).
# A run counts as finished only when the solver wrote "Computation time" into parameters.txt
# (all frames written is not enough: the solver can hang after the last one). Runs from
# before the run date was recorded have no such line; for them all frames is enough.
# fetch refuses to write into an existing local results/MyRun; --force replaces it.
# fetch without a name downloads every run on the cluster except the queued or running ones;
# runs that already exist locally are skipped, or replaced with --force.
# cancel also deletes the upload folder of a job that has not started yet, so its name is free again.
# clean deletes the upload folders in jobs/ that no queued or running job will use (leftovers,
# e.g. from a job that was killed while it was starting).
# delete-results deletes results/ and jobs/ (so the run names can be used again). It refuses
# while you have a job queued or running, lists the runs and asks you to type "yes";
# --yes skips the question. Your local results are kept.
#
# Configure with environment variables or by editing the defaults below.
# Best is an alias in ~/.ssh/config, then CLUSTER_HOST=that alias.

set -euo pipefail

CLUSTER_HOST="${CLUSTER_HOST:-gp3}"                  # ssh alias from ~/.ssh/config
# Relative to remote $HOME. A dedicated folder; it holds results/, jobs/ and the slurm-*.out logs.
REMOTE_DIR="${REMOTE_DIR:-Phase_field_slurm}"

# Git Bash's own ssh cannot use the Windows ssh-agent (where the key is loaded),
# so prefer Windows OpenSSH when it exists. No path mangling for native programs.
export MSYS_NO_PATHCONV=1
SSH=ssh
[ -x /c/Windows/System32/OpenSSH/ssh.exe ] && SSH=/c/Windows/System32/OpenSSH/ssh.exe

# ssh for commands that need no input. Without </dev/null, Windows ssh.exe started from
# Git Bash can wait for the terminal's input forever and the script hangs.
rssh() { "$SSH" "$CLUSTER_HOST" "$@" </dev/null; }

cd "$(git rev-parse --show-toplevel)"

# Run names end up in paths and in `rm -rf`, so allow only safe characters.
check_name() {
    case "$1" in
        ""|.|..|*[!A-Za-z0-9._-]*)
            echo "Invalid run name '$1': use only letters, digits, '.', '_' and '-'." >&2
            exit 1 ;;
    esac
}

cmd="${1:-}"; shift || true

case "$cmd" in
submit|test)
    # test is submit with a second archive (the baseline) and MODE=test for run.sbatch.
    mode=$cmd
    folder=""; delete=0; force=0; gpu=""; against=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --delete-calcs)
                [ "$mode" = submit ] || { echo "--delete-calcs is only for submit." >&2; exit 1; }
                delete=1 ;;
            --force)        force=1 ;;
            --gpu)
                [ $# -ge 2 ] || { echo "--gpu needs a GPU type, e.g. --gpu rtx5080" >&2; exit 1; }
                gpu="$2"; shift ;;
            --against)
                [ "$mode" = test ] || { echo "--against is only for test." >&2; exit 1; }
                [ $# -ge 2 ] || { echo "--against needs a commit, e.g. --against HEAD~1" >&2; exit 1; }
                against="$2"; shift ;;
            -*) echo "Unknown option: $1" >&2; exit 1 ;;
            *)  folder="$1" ;;
        esac
        shift
    done
    if [ -z "$folder" ]; then
        if [ "$mode" = test ]; then
            echo "Usage: cluster.sh test NAME [--against COMMIT] [--force] [--gpu TYPE]" >&2
        else
            echo "Usage: cluster.sh submit FOLDER [--delete-calcs] [--force] [--gpu TYPE]" >&2
        fi
        exit 1
    fi
    check_name "$folder"
    job_dir="jobs/$folder"
    if [ "$mode" = submit ]; then
        [ -f config/config.json ] || { echo "config/config.json does not exist, nothing to run." >&2; exit 1; }
    else
        [ -n "$(ls tests/configs/*.json 2>/dev/null)" ] || { echo "There are no tests/configs/*.json, nothing to test." >&2; exit 1; }
        against="${against:-HEAD}"
        base_commit=$(git rev-parse --verify --quiet "$against^{commit}") \
            || { echo "'$against' is not a commit." >&2; exit 1; }
        base_short=$(git rev-parse --short "$base_commit")
        # The test compares the working tree with the baseline. If they have the same code (e.g.
        # the changes were committed before testing and the baseline is HEAD), both versions are
        # the same and every case would come out IDENTICAL without saying anything.
        # config/config.json does not count: the test does not use it.
        if git diff --quiet "$base_commit" -- . ':(exclude)config/config.json' \
           && [ -z "$(git ls-files -o --exclude-standard)" ]; then
            echo "Your working tree has the same code as $against ($base_short), so the test would compare the code with itself." >&2
            echo "If you have already committed your changes, compare with the commit before them, e.g. --against HEAD~1." >&2
            exit 1
        fi
    fi

    # --gpu TYPE requests that card type (--gres=gpu:TYPE:1 overrides the "any GPU" line in
    # run.sbatch). Checked against the cluster before anything is uploaded.
    gres=()
    if [ -n "$gpu" ]; then
        case "$gpu" in
            *[!A-Za-z0-9_]*) echo "Invalid GPU type '$gpu'." >&2; exit 1 ;;
        esac
        types=$(rssh "sinfo -h -o %G" | tr ',' '\n' | sed -nE 's/^gpu:([^:]+):.*/\1/p' | sort -u)
        if ! printf '%s\n' "$types" | grep -qx "$gpu"; then
            echo "GPU type '$gpu' does not exist on $CLUSTER_HOST. Available types:" $types >&2
            exit 1
        fi
        gres=(--gres="gpu:$gpu:1")
    fi

    # A name is taken if its results exist (finished or running) or its upload waits in jobs/ (queued).
    taken=$(rssh "cd '$REMOTE_DIR' 2>/dev/null && { [ -e 'results/$folder' ] || [ -e '$job_dir' ]; } && echo yes || echo no")
    if [ "$taken" = "yes" ]; then
        if [ "$force" = "1" ]; then
            echo ">> '$folder' already exists on $CLUSTER_HOST: deleting old results and upload (--force)"
            rssh "cd '$REMOTE_DIR' && rm -rf 'results/$folder' '$job_dir'"
        else
            echo "Run name '$folder' is already used on $CLUSTER_HOST (finished or queued)." >&2
            echo "Pick another name, or use --force to delete the old run and start again." >&2
            exit 1
        fi
    fi

    echo ">> Uploading sources to $CLUSTER_HOST:$REMOTE_DIR/$job_dir/sources.tar.gz"
    # Tracked + untracked-but-not-ignored files (with config/config.json); results/ and build/
    # are gitignored. Files deleted locally but still in the index would break tar, so filter them.
    # The archive belongs to this run only: the job unpacks it into its own folder and builds
    # it from scratch, so nothing from other runs can get mixed in.
    git ls-files -co --exclude-standard -z \
        | while IFS= read -r -d '' f; do [ -e "$f" ] && printf '%s\0' "$f"; done \
        | tar --null -T - -czf - \
        | "$SSH" "$CLUSTER_HOST" "mkdir -p '$REMOTE_DIR/$job_dir' && cat > '$REMOTE_DIR/$job_dir/sources.tar.gz'"

    if [ "$mode" = test ]; then
        # The committed baseline, exactly as git stores it. core.autocrlf=false: old commits have
        # no .gitattributes, and CRLF line endings would break make and the scripts on Linux.
        echo ">> Uploading the baseline $against ($base_short) to $CLUSTER_HOST:$REMOTE_DIR/$job_dir/baseline.tar.gz"
        git -c core.autocrlf=false archive --format=tar.gz "$base_commit" \
            | "$SSH" "$CLUSTER_HOST" "cat > '$REMOTE_DIR/$job_dir/baseline.tar.gz'"
    fi

    # Local commit, marked "-dirty" when uncommitted changes (or untracked files) were uploaded too.
    # config/config.json does not count: it is the input of the run, not code, and every run
    # keeps it in info/setup/config.json (and in sources.tar.gz) anyway.
    commit=$(git rev-parse --short HEAD)
    [ -z "$(git status --porcelain -- . ':(exclude)config/config.json')" ] || commit="$commit-dirty"

    # The batch script comes from the uploaded archive too (sbatch reads it from stdin), so the
    # job script always matches the code. The job name is the run name, which `status` shows
    # and `cancel` and `clean` rely on. If sbatch fails, the upload is deleted so the name stays free.
    # A test takes minutes, so it gets a short time limit (overrides the 12 h in run.sbatch):
    # if something hangs, it does not hold the GPU for long.
    extra=()
    if [ "$mode" = test ]; then
        extra=(--time=01:00:00)
        echo ">> Submitting test job (commit $commit against the baseline $base_short, GPU ${gpu:-any})"
    else
        echo ">> Submitting job (commit $commit, GPU ${gpu:-any})"
    fi
    rssh "cd '$REMOTE_DIR' && { tar xzOf '$job_dir/sources.tar.gz' slurm/run.sbatch | sbatch --job-name='$folder' ${gres[*]:-} ${extra[*]:-} --export=ALL,FOLDER='$folder',DELETE_CALCS=$delete,GIT_COMMIT='$commit',MODE=$mode,BASELINE_COMMIT='${base_short:-}' || { rm -rf '$job_dir'; exit 1; }; }"
    ;;
status)
    # Queue plus progress of running jobs, read from the "Steps completed" lines that the
    # solver prints into slurm-<jobid>.out. The computation is finished only when the solver
    # has written "Computation time" into parameters.txt (see `list`).
    # The remote script is sent on stdin.
    "$SSH" "$CLUSTER_HOST" bash -s -- "$REMOTE_DIR" <<'REMOTE_SCRIPT'
cd "$1" 2>/dev/null || { echo "Folder $1 does not exist on the cluster yet (nothing was submitted)."; exit 0; }
# Seconds as "3 min 20 s".
duration() {
    if [ "$1" -ge 60 ]; then echo "$(( $1 / 60 )) min $(( $1 % 60 )) s"; else echo "$1 s"; fi
}
# Progress of the regression test in results/$2 (job $1), from the progress file that
# tests/regression.sh rewrites at every step. No estimate of the time left: the cases take very
# different times.
test_progress() {
    local phase a b c d e
    phase=""
    [ -f "results/$2/progress" ] && read -r phase a b c d e < "results/$2/progress"
    case "$phase" in
        build)   echo "job $1: regression test, compiling the $a version" ;;
        compare) echo "job $1: regression test, comparing the results" ;;
        run)
            # a = this run, b = number of runs, c = start of this run, d = case, e = version.
            echo "job $1: regression test, run $a/$b ($d, $e version), running for $(duration $(( $(date +%s) - c )))" ;;
        # Tests submitted before regression.sh wrote the progress file, or the job is still
        # unpacking the sources.
        *)       echo "job $1: regression test running (see slurm-$1.out)" ;;
    esac
}

squeue -u "$USER"
squeue -u "$USER" -h -t RUNNING -o '%i %j' | while read -r id name; do
    # A regression test runs many short cases, so it reports which run it is at.
    if [ -f "results/$name/info/setup/baseline.tar.gz" ]; then
        test_progress "$id" "$name"
        continue
    fi
    if grep -q '^Computation time:' "results/$name/info/setup/parameters.txt" 2>/dev/null; then
        echo "job $id: computation finished, now plotting"
        continue
    fi
    line=$(grep 'Steps completed' "slurm-$id.out" 2>/dev/null | tail -n 1)
    if [ -z "$line" ]; then
        echo "job $id: no progress yet (compiling or starting)"
        continue
    fi
    pct=$(echo "$line" | sed -E 's/.*=> *([0-9.]+)%.*/\1/')
    rem=$(echo "$line" | sed -E 's/.*Time remaining: *(.*[0-9]s) +Elapsed.*/\1/' | tr -s ' ')
    case "$pct" in
        100.00) echo "job $id: all frames written, but the solver has not finished (stuck? see slurm-$id.out)" ;;
        *)      echo "job $id: ${pct}% done, about ${rem} left" ;;
    esac
done
REMOTE_SCRIPT
    ;;
list)
    # One line per run in results/ or jobs/. The remote script is sent on stdin.
    "$SSH" "$CLUSTER_HOST" bash -s -- "$REMOTE_DIR" <<'REMOTE_SCRIPT'
cd "$1" 2>/dev/null || { echo "Folder $1 does not exist on the cluster yet (nothing was submitted)."; exit 0; }
shopt -s nullglob

# Your queued and running jobs, one per line: "name state jobid node" (the job name is the run name).
active=$(squeue -u "$USER" -h -o '%j %T %i %N') || exit 1

# Progress of running job $1 with parameters file $2, as in `status`: the percentage comes from
# the last "Steps completed" line in its log. The solver is done only when it has written
# "Computation time" into parameters.txt; 100% without it means it still runs after the last frame.
progress() {
    local line pct rem
    if [ -n "$2" ] && grep -q '^Computation time:' "$2" 2>/dev/null; then
        echo "computed, now plotting"; return
    fi
    line=$(grep 'Steps completed' "slurm-$1.out" 2>/dev/null | tail -n 1)
    if [ -z "$line" ]; then echo "compiling or starting"; return; fi
    pct=$(echo "$line" | sed -E 's/.*=> *([0-9.]+)%.*/\1/')
    rem=$(echo "$line" | sed -E 's/.*Time remaining: *(.*[0-9]s) +Elapsed.*/\1/' | tr -s ' ')
    case "$pct" in
        100.00) echo "all frames written, but the solver has not finished (stuck? see slurm-$1.out)" ;;
        *)      echo "${pct}% computed, about ${rem} left" ;;
    esac
}

# Newest log of a run: new jobs print "Run NAME, job ...", older ones only make's command line.
find_log() {
    local logs=(slurm-*.out)
    [ ${#logs[@]} -gt 0 ] || return 0
    grep -l -E -e "^Run $1, job " -e "\./results/$1 " "${logs[@]}" 2>/dev/null | tail -n 1
}

printf '%-22s %-9s %-9s %-16s %6s  %s\n' RUN STATE FRAMES "RUN DATE" SIZE DETAILS
names=$( { ls -1 results 2>/dev/null; ls -1 jobs 2>/dev/null; } | sort -u )
for name in $names; do
    r="results/$name"
    read -r jstate jid node <<< "$(printf '%s\n' "$active" | awk -v n="$name" '$1 == n { print $2, $3, $4; exit }')"

    # Results from before the info/ subfolders have parameters.txt and config.json directly in info/.
    params=""; cfg=""
    for f in "$r/info/setup/parameters.txt" "$r/info/parameters.txt"; do [ -f "$f" ] && { params=$f; break; }; done
    for f in "$r/info/setup/config.json" "$r/info/config.json"; do [ -f "$f" ] && { cfg=$f; break; }; done

    # Frames: the highest Result-NNNNN written so far / frame_num from the config.
    last=$(ls "$r/calculations" 2>/dev/null | sed -nE 's/^Result-([0-9]+)\..*/\1/p' | sort -n | tail -n 1)
    [ -n "$last" ] && last=$((10#$last))
    total=""
    [ -n "$cfg" ] && total=$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1]))["solver"]["frame_num"])' "$cfg" 2>/dev/null)
    frames="-"
    [ -n "$last$total" ] && frames="${last:-0}/${total:-?}"

    date="-"; comp_time=""; d=""
    if [ -n "$params" ]; then
        d=$(sed -n 's/^Run date: *//p' "$params" | cut -c1-16)
        [ -n "$d" ] && date=$d
        comp_time=$(sed -n 's/^Computation time: *//p' "$params" | sed 's/ *(.*//')
    fi
    # Since the binary records the run date, it also records the computation time when it ends,
    # so only that line means "finished". Older runs have neither; for them all frames have to do.
    complete=0
    if [ -n "$comp_time" ]; then
        complete=1
    elif [ -z "$d" ] && [ -n "$last" ] && [ "$last" = "$total" ]; then
        complete=1
    fi
    size="-"
    [ -d "$r" ] && size=$(du -sh "$r" | cut -f1)
    # A regression test (cluster.sh test) keeps the baseline next to the sources; it has no frames
    # of its own (they are in new/ and base/), and report.txt is written at its end.
    is_test=0
    [ -f "$r/info/setup/baseline.tar.gz" ] && is_test=1

    case "$jstate" in
        PENDING) state=queued;  details="job $jid" ;;
        RUNNING)
            state=running
            if [ "$is_test" = 1 ]; then details="job $jid on $node: regression test"
            else details="job $jid on $node: $(progress "$jid" "$params")"; fi ;;
        "")
            if [ "$is_test" = 1 ]; then
                if [ -f "$r/report.txt" ]; then
                    state=tested; details="$(tail -n 1 "$r/report.txt")"
                else
                    state=stopped
                    log=$(find_log "$name")
                    details="test did not finish${log:+; see $log}"
                fi
            elif [ ! -d "$r" ]; then
                # A job that failed before it moved its upload to the results left this behind.
                log=$(find_log "$name")
                state=leftover; details="upload in jobs/ without a job${log:+ (it failed at the start, see $log)}; run clean"
            elif [ "$complete" = 1 ]; then
                state=finished; details=${comp_time:+"computed in $comp_time"}
                # Animations are .mp4, or .gif without ffmpeg; old results have them directly in info/.
                if [ -z "$(find "$r/info" \( -name '*.mp4' -o -name '*.gif' \) -print -quit 2>/dev/null)" ]; then
                    details="${details:+$details; }no animations (plotting failed?)"
                fi
            else
                state=stopped
                log=$(find_log "$name")
                details="crashed, cancelled or time limit${log:+; see $log}"
                if [ -n "$last" ] && [ "$last" = "$total" ]; then
                    details="all frames written, but the solver did not finish: $details"
                fi
            fi ;;
        *) state=$(echo "$jstate" | tr '[:upper:]' '[:lower:]'); details="job $jid" ;;
    esac
    printf '%-22s %-9s %-9s %-16s %6s  %s\n' "$name" "$state" "$frames" "$date" "$size" "$details"
done
REMOTE_SCRIPT
    ;;
log)
    rssh "cd '$REMOTE_DIR' && tail -n 40 \"\$(ls -t slurm-*.out | head -1)\""
    ;;
gpus)
    # GPU types on the cluster and how many are free right now (from sinfo), plus the pending
    # jobs of all users that wait for each type (from squeue). The remote script is sent on stdin.
    "$SSH" "$CLUSTER_HOST" bash -s <<'REMOTE_SCRIPT'
{
    # Node lines: "node gpu:TYPE:COUNT(S:0) gpu:TYPE:USED(IDX:..) state"
    sinfo -N -h -p gp -O "NodeList:40,Gres:80,GresUsed:80,StateLong:40" | sed 's/^/node /'
    # Pending jobs: "pending gres/gpu:TYPE:1" or "pending gres/gpu:1" (any type)
    squeue -h -t PENDING -o "%b" | sed 's/^/pending /'
} | awk '
$1 == "node" && $3 ~ /^gpu:/ {
    split($3, g, ":"); type = g[2]; total = g[3] + 0     # "gpu:rtx3060:2(S:0)"   -> 2
    split($4, u, ":"); used = u[3] + 0                   # "gpu:rtx3060:1(IDX:0)" -> 1
    nodes[type] = nodes[type] (nodes[type] == "" ? "" : ",") $2
    tot[type] += total; usd[type] += used
    # GPUs on a drained/down node cannot get new jobs
    if ($5 ~ /drain|down|fail|maint|reserv|not_resp|unk/) { unav[type] += total - used; nodes[type] = nodes[type] "(" $5 ")" }
}
$1 == "pending" {
    n = split($2, p, ":")
    if (n == 3) wait[p[2]]++                             # gres/gpu:TYPE:1
    else if ($2 ~ /gpu/) wait_any++                      # gres/gpu:1
}
END {
    printf "%-9s %5s %7s %12s %8s  %s\n", "TYPE", "FREE", "IN USE", "UNAVAILABLE", "WAITING", "NODES"
    cmd = "sort"
    for (t in tot)
        printf "%-9s %5d %7d %12d %8d  %s\n", t, tot[t] - usd[t] - unav[t], usd[t], unav[t], wait[t], nodes[t] | cmd
    close(cmd)
    printf "\nPending jobs of all users that accept any GPU type: %d\n", wait_any
    print "UNAVAILABLE = idle GPUs on drained or down nodes; WAITING = pending jobs that asked for that type."
    print "Choose a type with: cluster.cmd submit NAME --gpu TYPE"
}'
REMOTE_SCRIPT
    ;;
fetch)
    folder=""; no_calcs=0; force=0
    for arg in "$@"; do
        case "$arg" in
            --no-calcs) no_calcs=1 ;;
            --force)    force=1 ;;
            -*) echo "Unknown option: $arg" >&2; exit 1 ;;
            *)  folder="$arg" ;;
        esac
    done

    # Downloads results/$1 from the cluster into the local results/; an existing local copy is
    # replaced first (the callers decide whether that is allowed).
    fetch_one() {
        local f=$1 excl=""
        if [ -e "results/$f" ]; then
            echo ">> Replacing existing local results/$f (--force)"
            rm -rf "results/$f"
        fi
        # Excluded when packing on the cluster, so the calculations are not downloaded at all.
        # (GNU tar applies --exclude only to the names after it.)
        [ "$no_calcs" = "1" ] && excl="--exclude='$f/calculations'"
        mkdir -p results
        rssh "cd '$REMOTE_DIR/results' && tar czf - $excl '$f'" | tar xzf - -C results
        echo ">> Saved to results/$f"
    }

    if [ -n "$folder" ]; then
        check_name "$folder"
        if [ -e "results/$folder" ] && [ "$force" != "1" ]; then
            echo "Local results/$folder already exists. Remove it or use --force to replace it." >&2
            exit 1
        fi
        fetch_one "$folder"
        exit 0
    fi

    # No name: every run on the cluster, except the ones that are still queued or running
    # (their results are incomplete) and, without --force, the ones you already have locally.
    runs=$(rssh "cd '$REMOTE_DIR/results' 2>/dev/null && ls -1A || true")
    if [ -z "$runs" ]; then
        echo "There are no results on $CLUSTER_HOST:$REMOTE_DIR, nothing to fetch."
        exit 0
    fi
    active=$(rssh "squeue -u \"\$USER\" -h -o %j")      # job name = run name
    fetched=0; skipped=0
    for f in $runs; do
        case "$f" in
            .|..|*[!A-Za-z0-9._-]*)
                echo ">> Skipping '$f': not a valid run name"; skipped=$((skipped + 1)); continue ;;
        esac
        if printf '%s\n' "$active" | grep -qxF -- "$f"; then
            echo ">> Skipping $f: still queued or running"; skipped=$((skipped + 1)); continue
        fi
        if [ -e "results/$f" ] && [ "$force" != "1" ]; then
            echo ">> Skipping $f: results/$f already exists locally (--force replaces it)"
            skipped=$((skipped + 1)); continue
        fi
        fetch_one "$f"
        fetched=$((fetched + 1))
    done
    echo ">> Fetched $fetched run(s), skipped $skipped"
    ;;
cancel)
    id="${1:?Usage: cluster.sh cancel JOBID}"
    case "$id" in *[!0-9]*) echo "Invalid job ID '$id'." >&2; exit 1 ;; esac
    # A job that has not started yet still has its upload in jobs/<name> (the job name is the
    # run name), which would keep the name taken, so delete it. A running job has already moved
    # its sources to the results and cleans up after itself. The remote script is sent on stdin.
    "$SSH" "$CLUSTER_HOST" bash -s -- "$REMOTE_DIR" "$id" <<'REMOTE_SCRIPT'
read -r state name <<< "$(squeue -h -j "$2" -o '%T %j' 2>/dev/null)"
scancel "$2" || exit 1
case "$name" in ""|.|..|*[!A-Za-z0-9._-]*) exit 0 ;; esac
if [ "$state" = PENDING ] && [ -d "$1/jobs/$name" ]; then
    rm -rf "$1/jobs/$name"
    echo ">> Deleted the upload jobs/$name of the cancelled job"
fi
REMOTE_SCRIPT
    ;;
clean)
    # Upload folders whose name is not the name of a queued or running job are leftovers.
    # The remote script is sent on stdin.
    "$SSH" "$CLUSTER_HOST" bash -s -- "$REMOTE_DIR" <<'REMOTE_SCRIPT'
cd "$1/jobs" 2>/dev/null || { echo "There is no jobs/ folder on the cluster, nothing to clean."; exit 0; }
active=$(squeue -u "$USER" -h -o %j) || exit 1
deleted=0
for d in *; do
    [ -d "$d" ] || continue
    printf '%s\n' "$active" | grep -qxF -- "$d" && continue
    rm -rf -- "$d"
    echo ">> Deleted leftover upload jobs/$d"
    deleted=$((deleted + 1))
done
echo ">> $deleted leftover upload(s) deleted"
REMOTE_SCRIPT
    ;;
delete-results)
    yes=0
    for arg in "$@"; do
        case "$arg" in
            --yes) yes=1 ;;
            *) echo "Unknown option: $arg" >&2; exit 1 ;;
        esac
    done
    # A running job would lose its output and a queued job its sources, so there is no override.
    active=$(rssh "squeue -u \"\$USER\" -h | wc -l" | tr -d '[:space:]')
    if [ "$active" != "0" ]; then
        echo "You have $active job(s) queued or running on $CLUSTER_HOST (see 'status')." >&2
        echo "Wait for them to finish (or cancel them) before deleting the results." >&2
        exit 1
    fi

    runs=$(rssh "cd '$REMOTE_DIR/results' 2>/dev/null && ls -1A")
    if [ -z "$runs" ]; then
        echo "There are no results on $CLUSTER_HOST:$REMOTE_DIR, nothing to delete."
        exit 0
    fi
    size=$(rssh "du -sh '$REMOTE_DIR/results' | cut -f1")
    echo "These $(echo "$runs" | wc -l | tr -d ' ') run(s) in $CLUSTER_HOST:$REMOTE_DIR/results ($size) will be deleted:"
    echo "$runs" | sed 's/^/    /'
    if [ "$yes" != "1" ]; then
        answer=""
        read -r -p "Type 'yes' to delete them: " answer || true
        if [ "$answer" != "yes" ]; then
            echo "Aborted, nothing was deleted."
            exit 1
        fi
    fi
    rssh "cd '$REMOTE_DIR' && rm -rf results jobs"
    echo ">> Deleted all results and uploads in $CLUSTER_HOST:$REMOTE_DIR"
    ;;
*)
    sed -n '2,37p' "$0"
    exit 1
    ;;
esac
