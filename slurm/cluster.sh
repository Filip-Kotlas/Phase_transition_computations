#!/bin/bash
# Run from your own computer (Git Bash / WSL / Linux), from anywhere in the repo.
#
#   slurm/cluster.sh submit MyRun [--delete-calcs] [--force]   upload code, sbatch the job
#   slurm/cluster.sh status                                    show your queue
#   slurm/cluster.sh log                                       tail the newest slurm-*.out
#   slurm/cluster.sh fetch MyRun [--no-calcs] [--force]        download results/MyRun
#   slurm/cluster.sh cancel JOBID
#
# submit refuses a name that is already used on the cluster (finished or still queued);
# --force deletes the old results and config of that name first. It does NOT stop a job
# that is still running under that name, so check `status` first.
# fetch refuses to write into an existing local results/MyRun; --force replaces it.
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
    folder=""; delete=0; force=0
    for arg in "$@"; do
        case "$arg" in
            --delete-calcs) delete=1 ;;
            --force)        force=1 ;;
            -*) echo "Unknown option: $arg" >&2; exit 1 ;;
            *)  folder="$arg" ;;
        esac
    done
    [ -n "$folder" ] || { echo "Usage: cluster.sh submit FOLDER [--delete-calcs] [--force]" >&2; exit 1; }
    check_name "$folder"
    remote_cfg="config/queued/$folder.json"

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

    echo ">> Submitting job (commit $commit)"
    rssh "cd '$REMOTE_DIR' && sbatch --export=ALL,FOLDER='$folder',CONFIG='$remote_cfg',DELETE_CALCS=$delete,GIT_COMMIT='$commit' slurm/run.sbatch"
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
*)
    sed -n '2,15p' "$0"
    exit 1
    ;;
esac
