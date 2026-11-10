#!/bin/bash
# Run from your own computer (Git Bash / WSL / Linux), from anywhere in the repo.
#
#   slurm/cluster.sh submit MyRun [--delete-calcs] [--force] [--gpu TYPE]
#                                                              upload code, sbatch the job
#                                                              (--gpu rtx5080: only that card type)
#   slurm/cluster.sh status                                    show your queue
#   slurm/cluster.sh log                                       tail the newest slurm-*.out
#   slurm/cluster.sh gpus                                      GPU types and how many are free now
#   slurm/cluster.sh fetch MyRun [--no-calcs] [--force]        download results/MyRun
#   slurm/cluster.sh cancel JOBID
#   slurm/cluster.sh clean [--force]                           make clean on the cluster (deletes build/)
#   slurm/cluster.sh delete-results [--yes]                    delete ALL results on the cluster
#
# submit refuses a name that is already used on the cluster (finished or still queued);
# --force deletes the old results and config of that name first. It does NOT stop a job
# that is still running under that name, so check `status` first.
# fetch refuses to write into an existing local results/MyRun; --force replaces it.
# clean refuses while you have a job queued or running (it could be compiling right now);
# --force cleans anyway. Needed after changing compiler flags in the Makefile.
# delete-results deletes results/ and the per-run configs in config/queued/ (so the run
# names can be used again). It refuses while you have a job queued or running, lists the
# runs and asks you to type "yes"; --yes skips the question. Your local results are kept.
#
# Configure with environment variables or by editing the defaults below.
# Best is an alias in ~/.ssh/config, then CLUSTER_HOST=that alias.

set -euo pipefail

CLUSTER_HOST="${CLUSTER_HOST:-gp3}"                  # ssh alias from ~/.ssh/config
# Relative to remote $HOME. A dedicated folder: submit replaces src/, include/ and python/ in it.
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
submit)
    folder=""; delete=0; force=0; gpu=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --delete-calcs) delete=1 ;;
            --force)        force=1 ;;
            --gpu)
                [ $# -ge 2 ] || { echo "--gpu needs a GPU type, e.g. --gpu rtx5080" >&2; exit 1; }
                gpu="$2"; shift ;;
            -*) echo "Unknown option: $1" >&2; exit 1 ;;
            *)  folder="$1" ;;
        esac
        shift
    done
    [ -n "$folder" ] || { echo "Usage: cluster.sh submit FOLDER [--delete-calcs] [--force] [--gpu TYPE]" >&2; exit 1; }
    check_name "$folder"
    remote_cfg="config/queued/$folder.json"

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

    # A name is taken if its results exist or its config was already queued.
    taken=$(rssh "cd '$REMOTE_DIR' 2>/dev/null && { [ -e 'results/$folder' ] || [ -e '$remote_cfg' ]; } && echo yes || echo no")
    if [ "$taken" = "yes" ]; then
        if [ "$force" = "1" ]; then
            echo ">> '$folder' already exists on $CLUSTER_HOST: deleting old results and config (--force)"
            rssh "cd '$REMOTE_DIR' && rm -rf 'results/$folder' '$remote_cfg'"
        else
            echo "Run name '$folder' is already used on $CLUSTER_HOST (finished or queued)." >&2
            echo "Pick another name, or use --force to delete the old run and start again." >&2
            exit 1
        fi
    fi

    echo ">> Uploading sources to $CLUSTER_HOST:$REMOTE_DIR"
    # Tracked + untracked-but-not-ignored files; results/ and build/ are gitignored.
    # Files deleted locally but still in the index would break tar, so filter them.
    # src/, include/ and python/ are wiped remotely first, so a file you deleted or renamed
    # locally cannot linger and get compiled (the Makefile builds every src/*.cu|cpp|c).
    # tar keeps modification times, so unchanged files do not trigger a recompile.
    git ls-files -co --exclude-standard -z \
        | while IFS= read -r -d '' f; do [ -e "$f" ] && printf '%s\0' "$f"; done \
        | tar --null -T - -czf - \
        | "$SSH" "$CLUSTER_HOST" "mkdir -p '$REMOTE_DIR' && rm -rf '$REMOTE_DIR/src' '$REMOTE_DIR/include' '$REMOTE_DIR/python' && tar xzf - -C '$REMOTE_DIR'"

    # Per-run copy of the config, so later submits cannot change a queued job.
    echo ">> Uploading config/config.json as $remote_cfg"
    "$SSH" "$CLUSTER_HOST" "mkdir -p '$REMOTE_DIR/config/queued' && cat > '$REMOTE_DIR/$remote_cfg'" < config/config.json

    # Local commit, marked "-dirty" when uncommitted changes were uploaded too.
    commit=$(git rev-parse --short HEAD)
    [ -z "$(git status --porcelain)" ] || commit="$commit-dirty"

    echo ">> Submitting job (commit $commit, GPU ${gpu:-any})"
    rssh "cd '$REMOTE_DIR' && sbatch ${gres[*]:-} --export=ALL,FOLDER='$folder',CONFIG='$remote_cfg',DELETE_CALCS=$delete,GIT_COMMIT='$commit' slurm/run.sbatch"
    ;;
status)
    # Queue plus progress of running jobs, read from the "Steps completed" lines that the
    # solver prints into slurm-<jobid>.out. The remote script is sent on stdin.
    "$SSH" "$CLUSTER_HOST" bash -s -- "$REMOTE_DIR" <<'REMOTE_SCRIPT'
cd "$1" 2>/dev/null || { echo "Folder $1 does not exist on the cluster yet (nothing was submitted)."; exit 0; }
squeue -u "$USER"
for id in $(squeue -u "$USER" -h -t RUNNING -o %i); do
    line=$(grep 'Steps completed' "slurm-$id.out" 2>/dev/null | tail -n 1)
    if [ -z "$line" ]; then
        echo "job $id: no progress yet (compiling or starting)"
        continue
    fi
    pct=$(echo "$line" | sed -E 's/.*=> *([0-9.]+)%.*/\1/')
    rem=$(echo "$line" | sed -E 's/.*Time remaining: *(.*[0-9]s) +Elapsed.*/\1/' | tr -s ' ')
    case "$pct" in
        100.00) echo "job $id: computation finished, now plotting" ;;
        *)      echo "job $id: ${pct}% done, about ${rem} left" ;;
    esac
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
    [ -n "$folder" ] || { echo "Usage: cluster.sh fetch FOLDER [--no-calcs] [--force]" >&2; exit 1; }
    check_name "$folder"

    if [ -e "results/$folder" ]; then
        if [ "$force" = "1" ]; then
            echo ">> Replacing existing local results/$folder (--force)"
            rm -rf "results/$folder"
        else
            echo "Local results/$folder already exists. Remove it or use --force to replace it." >&2
            exit 1
        fi
    fi

    excl=()
    [ "$no_calcs" = "1" ] && excl=(--exclude="$folder/calculations")
    mkdir -p results
    rssh "cd '$REMOTE_DIR/results' && tar czf - '$folder'" | tar xzf - -C results "${excl[@]}"
    echo ">> Saved to results/$folder"
    ;;
cancel)
    rssh "scancel '${1:?Usage: cluster.sh cancel JOBID}'"
    ;;
clean)
    force=0
    for arg in "$@"; do
        case "$arg" in
            --force) force=1 ;;
            *) echo "Unknown option: $arg" >&2; exit 1 ;;
        esac
    done
    # A job that is compiling would lose its objects halfway, so refuse while any job is active.
    active=$(rssh "squeue -u \"\$USER\" -h | wc -l" | tr -d '[:space:]')
    if [ "$active" != "0" ] && [ "$force" != "1" ]; then
        echo "You have $active job(s) queued or running on $CLUSTER_HOST (see 'status')." >&2
        echo "Wait for them to finish, or use --force to clean anyway." >&2
        exit 1
    fi
    rssh "[ -d '$REMOTE_DIR' ] || { echo 'Folder $REMOTE_DIR does not exist on the cluster yet, nothing to clean.'; exit 0; }; cd '$REMOTE_DIR' && make clean"
    echo ">> Cleaned build/ in $CLUSTER_HOST:$REMOTE_DIR"
    ;;
delete-results)
    yes=0
    for arg in "$@"; do
        case "$arg" in
            --yes) yes=1 ;;
            *) echo "Unknown option: $arg" >&2; exit 1 ;;
        esac
    done
    # A running job would lose its output and a queued job its config, so there is no override.
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
    rssh "cd '$REMOTE_DIR' && rm -rf results && rm -f config/queued/*.json"
    echo ">> Deleted all results and queued configs in $CLUSTER_HOST:$REMOTE_DIR"
    ;;
*)
    sed -n '2,25p' "$0"
    exit 1
    ;;
esac
