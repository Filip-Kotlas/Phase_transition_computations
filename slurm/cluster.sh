#!/bin/bash
# Run from your own computer (Git Bash / WSL / Linux), from anywhere in the repo.
#
#   slurm/cluster.sh submit MyRun [--config FILE] [--delete-calcs] [--force] [--gpu TYPE] [--mem SIZE]
#                                                              upload code, sbatch the job
#                                                              (--config: run FILE instead of
#                                                              config/config.json;
#                                                              --gpu rtx5080: only that card type;
#                                                              --mem 4G: host RAM instead of the
#                                                              default in run.sbatch; also for
#                                                              test and bench)
#   slurm/cluster.sh submit-many [FOLDER|FILE...] [--prefix P] [--dry-run] [--delete-calcs] [--force] [--gpu TYPE] [--mem SIZE]
#                                                              one submit per config: every *.json
#                                                              in FOLDER (default
#                                                              comp_configurations/) or the FILEs;
#                                                              run name = P + file name without
#                                                              .json (--dry-run: only the checks)
#   slurm/cluster.sh test MyTest [--against COMMIT] [--new COMMIT] [--force] [--gpu TYPE]
#                                                              regression test: your working tree
#                                                              (or --new COMMIT) against COMMIT
#                                                              (default HEAD)
#   slurm/cluster.sh bench MyBench [--against COMMIT] [--new COMMIT] [--repeat N] [--long] [--force] [--gpu TYPE]
#                                                              benchmark: computation times of your
#                                                              working tree (or --new COMMIT)
#                                                              against COMMIT (--long: one demanding
#                                                              case, about 10 min per run)
#   slurm/cluster.sh status                                   show your queue
#   slurm/cluster.sh list                                      all runs on the cluster and their state
#   slurm/cluster.sh log                                       tail the newest slurm-*.out
#   slurm/cluster.sh gpus [--mem SIZE]                         before submitting: free and usable GPUs
#                                                              per type, what blocks the free ones,
#                                                              your jobs, and whether a job starts now
#   slurm/cluster.sh nodes [--mem SIZE]                        every node: free GPUs, CPUs and RAM,
#                                                              its running jobs, and whether our
#                                                              job would fit there now
#   slurm/cluster.sh fetch [MyRun] [--no-calcs] [--force]      download results/MyRun, or all runs
#   slurm/cluster.sh cancel JOBID
#   slurm/cluster.sh clean                                     delete leftover upload folders jobs/*
#   slurm/cluster.sh delete-results [--yes]                    delete ALL results on the cluster
#
# Every run gets its own copy of the code: submit uploads the sources (including
# config/config.json) as jobs/MyRun/sources.tar.gz, and the job compiles and runs that copy
# on the node's local disk. A later submit therefore never changes the code of a queued job.
# The job keeps the archive as results/MyRun/info/setup/sources.tar.gz and deletes jobs/MyRun.
# submit --config FILE puts FILE into the archive as config/config.json instead of your
# config/config.json, which stays untouched (a relative FILE is relative to the folder you
# started the script from). submit-many does that for every config it is given, each as its own
# run and job, so status, list, log, cancel and fetch see ordinary runs. It checks everything
# first (the files are JSON, the names are valid, distinct and free on the cluster, the GPU type
# exists) and submits nothing if a check fails; --force deletes the taken names first, as with
# submit. Under PowerShell *.json is not expanded, so give a folder rather than a pattern.
# submit refuses a name that is already used on the cluster (finished or still queued);
# --force deletes the old results and upload of that name first. It does NOT stop a job
# that is still running under that name, so check `status` first.
# test uploads the working tree like submit, plus the committed version COMMIT (git archive) as
# the baseline. The job builds both, runs every tests/configs/*.json with both on the same GPU
# and compares them (tests/regression.sh, python/compare_runs.py); `log` shows the result and
# results/MyTest/report.txt keeps it. So test your changes BEFORE committing them: once they are
# committed, HEAD is the new code itself, and test refuses to compare the code with itself
# (use --against HEAD~1, or the commit before your changes, then).
# bench uploads the same two versions as test, runs every tests/benchmark_configs/*.json N times
# (default 3) with each of them on the same GPU and compares the computation times
# (tests/benchmark.sh, python/compare_times.py); `log` shows the result and
# results/MyBench/report.txt keeps it. Unlike test it also runs when both versions have the same
# code: that measures the noise of the timings.
# bench --long runs tests/benchmark_configs_long/*.json instead: one demanding case (400x400 with a
# thin interface) of about 10 minutes per run on an RTX 2070 SUPER, once with each version (unless
# --repeat N). It shows whether a speedup holds for a long run.
# --new COMMIT (test and bench) puts that commit in place of your working tree as the version
# under test, e.g. `bench Old --against 5e9bccb~1 --new 5e9bccb` for an old change. Your working
# tree is still uploaded: the job script, the test scripts, the configs and the comparison come
# from it, so they also work for commits from before they existed (the commits still have to
# accept `Phase-field RESULT_FOLDER CONFIG` and read today's configs).
# list shows each run in results/ and jobs/ with its state: queued, running (with progress),
# finished (computation time; warns when the animations are missing), failed (the solver stopped
# because of NaN or inf, got stuck or could not write a frame, and wrote the reason into
# parameters.txt) or stopped (no job any
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

# Config files on the command line are relative to the folder the script was started from.
start_dir="$(pwd)"
cd "$(git rev-parse --show-toplevel)"

# Run names end up in paths and in `rm -rf`, so allow only safe characters.
check_name() {
    case "$1" in
        ""|.|..|*[!A-Za-z0-9._-]*)
            echo "Invalid run name '$1': use only letters, digits, '.', '_' and '-'." >&2
            exit 1 ;;
    esac
}

# A Slurm memory size: a whole number with an optional unit K, M (the default), G or T.
check_mem() {
    [[ "$1" =~ ^[1-9][0-9]*[KMGT]?$ ]] \
        || { echo "Invalid memory size '$1': use e.g. 4G or 3000M." >&2; exit 1; }
}

# The same size in MB (to compare it with what sinfo reports).
mem_to_mb() {
    local n=${1%[KMGT]}
    case "$1" in
        *K) echo $(( n / 1024 )) ;;
        *G) echo $(( n * 1024 )) ;;
        *T) echo $(( n * 1024 * 1024 )) ;;
        *)  echo "$n" ;;
    esac
}

# A path from the command line as this bash can open it: backslashes (PowerShell) become
# slashes, and a relative path is taken relative to the folder the script was started from.
local_path() {
    local p="${1//\\//}"
    case "$p" in
        /*|[A-Za-z]:/*) printf '%s\n' "$p" ;;
        *)              printf '%s\n' "$start_dir/$p" ;;
    esac
}

# A config file to upload must exist and be valid JSON. The solver checks the rest, but only on
# the node; a broken file is cheaper to catch here. The file goes to Python on stdin, so the
# Windows Python needs no path conversion. Without Python the JSON check is skipped.
check_config_file() {
    [ -f "$1" ] || { echo "Config file $1 does not exist." >&2; exit 1; }
    command -v python3 >/dev/null 2>&1 || return 0
    local err
    if ! err=$(python3 -c 'import json, sys; json.load(sys.stdin)' < "$1" 2>&1); then
        echo "Config file $1 is not valid JSON: $(printf '%s\n' "$err" | tail -n 1)" >&2
        exit 1
    fi
}

# --gpu TYPE: checks that the card type exists on the cluster and sets gres for sbatch
# (--gres=gpu:TYPE:1 overrides the "any GPU" line in run.sbatch). Empty TYPE: any GPU.
set_gres() {
    gres=()
    [ -n "$1" ] || return 0
    case "$1" in
        *[!A-Za-z0-9_]*) echo "Invalid GPU type '$1'." >&2; exit 1 ;;
    esac
    local types
    types=$(rssh "sinfo -h -o %G" | tr ',' '\n' | sed -nE 's/^gpu:([^:]+):.*/\1/p' | sort -u)
    if ! printf '%s\n' "$types" | grep -qx "$1"; then
        echo "GPU type '$1' does not exist on $CLUSTER_HOST. Available types:" $types >&2
        exit 1
    fi
    gres=(--gres="gpu:$1:1")
}

# The sources of a run as a tar.gz on stdout: tracked + untracked-but-not-ignored files (with
# config/config.json); results/ and build/ are gitignored. Files deleted locally but still in
# the index would break tar, so they are filtered out. With a config file as the argument, that
# file goes into the archive as config/config.json instead (from a staging folder, so the local
# config/config.json is never touched).
pack_sources() {
    local config="${1:-}" stage=""
    local extra=()
    if [ -n "$config" ]; then
        stage=$(mktemp -d)
        mkdir -p "$stage/config"
        cp "$config" "$stage/config/config.json"
        extra=(-C "$stage" config/config.json)
    fi
    git ls-files -co --exclude-standard -z \
        | while IFS= read -r -d '' f; do
              [ -e "$f" ] || continue
              [ -n "$config" ] && [ "$f" = config/config.json ] && continue
              printf '%s\0' "$f"
          done \
        | tar --null -T - ${extra[@]+"${extra[@]}"} -czf -
    [ -z "$stage" ] || rm -rf "$stage"
}

# Uploads the sources of run $1 (with the config file $2, if given) to jobs/$1/sources.tar.gz.
upload_sources() {
    echo ">> Uploading sources to $CLUSTER_HOST:$REMOTE_DIR/jobs/$1/sources.tar.gz${2:+ (config $2)}"
    pack_sources "${2:-}" \
        | "$SSH" "$CLUSTER_HOST" "mkdir -p '$REMOTE_DIR/jobs/$1' && cat > '$REMOTE_DIR/jobs/$1/sources.tar.gz'"
}

# Local commit, marked "-dirty" when uncommitted changes (or untracked files) are uploaded too.
# config/config.json does not count: it is the input of the run, not code, and every run
# keeps it in info/setup/config.json (and in sources.tar.gz) anyway.
local_commit() {
    local c
    c=$(git rev-parse --short HEAD)
    [ -z "$(git status --porcelain -- . ':(exclude)config/config.json')" ] || c="$c-dirty"
    printf '%s\n' "$c"
}

# sbatch for run $1, whose archives are already uploaded. Reads the settings from the caller:
# gres, extra, mode, delete, commit, base_short, new_short, repeat, bench_set.
# The batch script comes from the uploaded archive too (sbatch reads it from stdin), so the
# job script always matches the code. The job name is the run name, which `status` shows
# and `cancel` and `clean` rely on. If sbatch fails, the upload is deleted so the name stays free.
# Prints sbatch's "Submitted batch job N".
# run.sbatch calls the mode of a plain submit "run" (and rejects "submit": from 6db179f until
# this mapping, every plain submit failed at the start of the job with "Unknown MODE 'submit'").
sbatch_run() {
    local job_mode=$mode
    [ "$job_mode" = submit ] && job_mode=run
    rssh "cd '$REMOTE_DIR' && { tar xzOf 'jobs/$1/sources.tar.gz' slurm/run.sbatch | sbatch --job-name='$1' ${gres[*]:-} ${extra[*]:-} --export=ALL,FOLDER='$1',DELETE_CALCS=$delete,GIT_COMMIT='$commit',MODE=$job_mode,BASELINE_COMMIT='${base_short:-}',NEW_COMMIT='${new_short:-}',REPEAT=${repeat:-},BENCH_SET='${bench_set:-}' || { rm -rf 'jobs/$1'; exit 1; }; }"
}

# For gpus and nodes: reads [--mem SIZE] and sets need_cpu and mem to what a job of ours asks
# for (run.sbatch, or --mem SIZE as with submit), and prints it.
our_job_request() {
    mem=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --mem)
                [ $# -ge 2 ] || { echo "--mem needs a size, e.g. --mem 4G" >&2; exit 1; }
                check_mem "$2"; mem="$2"; shift ;;
            *) echo "Usage: cluster.sh $cmd [--mem SIZE]" >&2; exit 1 ;;
        esac
        shift
    done
    need_cpu=$(sed -nE 's/^#SBATCH --cpus-per-task=([0-9]+).*/\1/p' slurm/run.sbatch)
    need_cpu=${need_cpu:-1}
    [ -n "$mem" ] || mem=$(sed -nE 's/^#SBATCH --mem=([0-9]+[KMGT]?).*/\1/p' slurm/run.sbatch)
    echo "Our job: 1 GPU, $need_cpu CPU(s), $mem RAM (change the RAM with --mem SIZE)"
    echo
}

cmd="${1:-}"; shift || true

case "$cmd" in
submit|test|bench)
    # test and bench are submit with a second archive (the baseline) and MODE=test or MODE=bench
    # for run.sbatch.
    mode=$cmd
    folder=""; delete=0; force=0; gpu=""; mem=""; against=""; new=""; repeat=""; long=0; config=""
    while [ $# -gt 0 ]; do
        case "$1" in
            --delete-calcs)
                [ "$mode" = submit ] || { echo "--delete-calcs is only for submit." >&2; exit 1; }
                delete=1 ;;
            --config)
                [ "$mode" = submit ] || { echo "--config is only for submit (test and bench use their own configs)." >&2; exit 1; }
                [ $# -ge 2 ] || { echo "--config needs a config file, e.g. --config comp_configurations/a.json" >&2; exit 1; }
                config=$(local_path "$2"); shift ;;
            --force)        force=1 ;;
            --gpu)
                [ $# -ge 2 ] || { echo "--gpu needs a GPU type, e.g. --gpu rtx5080" >&2; exit 1; }
                gpu="$2"; shift ;;
            --mem)
                [ $# -ge 2 ] || { echo "--mem needs a size, e.g. --mem 4G" >&2; exit 1; }
                check_mem "$2"
                mem="$2"; shift ;;
            --against)
                [ "$mode" != submit ] || { echo "--against is only for test and bench." >&2; exit 1; }
                [ $# -ge 2 ] || { echo "--against needs a commit, e.g. --against HEAD~1" >&2; exit 1; }
                against="$2"; shift ;;
            --new)
                [ "$mode" != submit ] || { echo "--new is only for test and bench." >&2; exit 1; }
                [ $# -ge 2 ] || { echo "--new needs a commit, e.g. --new HEAD" >&2; exit 1; }
                new="$2"; shift ;;
            --repeat)
                [ "$mode" = bench ] || { echo "--repeat is only for bench." >&2; exit 1; }
                [ $# -ge 2 ] || { echo "--repeat needs a number of runs, e.g. --repeat 5" >&2; exit 1; }
                case "$2" in ""|*[!0-9]*|0) echo "--repeat needs a positive number, not '$2'." >&2; exit 1 ;; esac
                repeat="$2"; shift ;;
            --long)
                [ "$mode" = bench ] || { echo "--long is only for bench." >&2; exit 1; }
                long=1 ;;
            -*) echo "Unknown option: $1" >&2; exit 1 ;;
            *)  folder="$1" ;;
        esac
        shift
    done
    if [ -z "$folder" ]; then
        case "$mode" in
            test)  echo "Usage: cluster.sh test NAME [--against COMMIT] [--new COMMIT] [--force] [--gpu TYPE] [--mem SIZE]" >&2 ;;
            bench) echo "Usage: cluster.sh bench NAME [--against COMMIT] [--new COMMIT] [--repeat N] [--long] [--force] [--gpu TYPE] [--mem SIZE]" >&2 ;;
            *)     echo "Usage: cluster.sh submit FOLDER [--config FILE] [--delete-calcs] [--force] [--gpu TYPE] [--mem SIZE]" >&2 ;;
        esac
        exit 1
    fi
    check_name "$folder"
    job_dir="jobs/$folder"
    # The long benchmark runs one case of about 10 minutes, once with each version by default.
    bench_set=""; [ "$long" = 1 ] && bench_set=long
    [ -n "$repeat" ] || { [ "$long" = 1 ] && repeat=1 || repeat=3; }
    if [ "$mode" = submit ]; then
        if [ -n "$config" ]; then
            check_config_file "$config"
        else
            [ -f config/config.json ] || { echo "config/config.json does not exist, nothing to run." >&2; exit 1; }
        fi
    else
        configs_dir=tests/configs
        [ "$mode" = bench ] && configs_dir=tests/benchmark_configs${bench_set:+_$bench_set}
        [ -n "$(ls $configs_dir/*.json 2>/dev/null)" ] || { echo "There are no $configs_dir/*.json, nothing to run." >&2; exit 1; }
        against="${against:-HEAD}"
        base_commit=$(git rev-parse --verify --quiet "$against^{commit}") \
            || { echo "'$against' is not a commit." >&2; exit 1; }
        base_short=$(git rev-parse --short "$base_commit")
        # The version under test: the working tree, or the commit given with --new.
        new_commit=""; new_short=""
        if [ -n "$new" ]; then
            new_commit=$(git rev-parse --verify --quiet "$new^{commit}") \
                || { echo "'$new' is not a commit." >&2; exit 1; }
            new_short=$(git rev-parse --short "$new_commit")
        fi
        # The test compares the version under test with the baseline. If they have the same code
        # (e.g. the changes were committed before testing and the baseline is HEAD), both versions
        # are the same and every case would come out IDENTICAL without saying anything.
        # config/config.json does not count: the test does not use it.
        # A benchmark of the same code is still useful: it shows the noise of the timings.
        if [ -n "$new_commit" ]; then
            git diff --quiet "$base_commit" "$new_commit" -- . ':(exclude)config/config.json' && same=1 || same=0
            tested="The commit $new ($new_short)"
        else
            git diff --quiet "$base_commit" -- . ':(exclude)config/config.json' \
                && [ -z "$(git ls-files -o --exclude-standard)" ] && same=1 || same=0
            tested="Your working tree"
        fi
        if [ "$same" = 1 ]; then
            if [ "$mode" = test ]; then
                echo "$tested has the same code as $against ($base_short), so the test would compare the code with itself." >&2
                [ -z "$new_commit" ] && echo "If you have already committed your changes, compare with the commit before them, e.g. --against HEAD~1." >&2
                exit 1
            fi
            echo ">> Note: $tested has the same code as $against ($base_short), so this benchmark measures only the noise of the timings."
        fi
    fi

    # Checked against the cluster before anything is uploaded.
    set_gres "$gpu"

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

    # The archive belongs to this run only: the job unpacks it into its own folder and builds
    # it from scratch, so nothing from other runs can get mixed in.
    upload_sources "$folder" "$config"

    if [ "$mode" != submit ]; then
        # The committed baseline, exactly as git stores it. core.autocrlf=false: old commits have
        # no .gitattributes, and CRLF line endings would break make and the scripts on Linux.
        echo ">> Uploading the baseline $against ($base_short) to $CLUSTER_HOST:$REMOTE_DIR/$job_dir/baseline.tar.gz"
        git -c core.autocrlf=false archive --format=tar.gz "$base_commit" \
            | "$SSH" "$CLUSTER_HOST" "cat > '$REMOTE_DIR/$job_dir/baseline.tar.gz'"
        # With --new, the version under test too; the working tree uploaded above then only
        # supplies the job script, the test scripts, the configs and the comparison.
        if [ -n "${new_commit:-}" ]; then
            echo ">> Uploading the version under test $new ($new_short) to $CLUSTER_HOST:$REMOTE_DIR/$job_dir/new.tar.gz"
            git -c core.autocrlf=false archive --format=tar.gz "$new_commit" \
                | "$SSH" "$CLUSTER_HOST" "cat > '$REMOTE_DIR/$job_dir/new.tar.gz'"
        fi
    fi

    # With --new, the version under test is exactly that commit.
    commit=$(local_commit)
    [ -n "${new_short:-}" ] && commit=$new_short

    # A test takes minutes, so it gets a short time limit (overrides the 12 h in run.sbatch):
    # if something hangs, it does not hold the GPU for long. A benchmark runs every case
    # 2 x N times for about half a minute each; 1 h for N = 3, plus 15 min for every further run.
    # The long benchmark runs its case 2 x N times for about 10 minutes: 30 min per repetition,
    # plus 30 min for the builds and as a margin.
    extra=()
    case "$mode" in
        test)
            extra=(--time=01:00:00)
            echo ">> Submitting test job (commit $commit against the baseline $base_short, GPU ${gpu:-any})" ;;
        bench)
            if [ "$long" = 1 ]; then minutes=$(( 30 + 30 * repeat ))
            else                     minutes=$(( 60 + 15 * (repeat > 3 ? repeat - 3 : 0) )); fi
            extra=(--time=$minutes)
            echo ">> Submitting ${bench_set:+$bench_set }benchmark job (commit $commit against the baseline $base_short, $repeat run(s) of each, GPU ${gpu:-any})" ;;
        *)
            echo ">> Submitting job (commit $commit, GPU ${gpu:-any})" ;;
    esac
    # --mem SIZE overrides the host RAM in run.sbatch. Slurm starts a job only on a node where
    # that much RAM is not yet reserved by other jobs, so a smaller request can start sooner
    # (see `nodes`).
    if [ -n "$mem" ]; then
        extra+=(--mem="$mem")
        echo ">> Host RAM: $mem"
    fi
    sbatch_run "$folder"
    ;;
submit-many)
    # One ordinary submit per config file (as with submit --config), named after the file.
    # Everything is checked before the first upload, so a mistake submits nothing.
    prefix=""; delete=0; force=0; gpu=""; mem=""; dry=0; inputs=()
    while [ $# -gt 0 ]; do
        case "$1" in
            --prefix)
                [ $# -ge 2 ] || { echo "--prefix needs a prefix, e.g. --prefix Conv_" >&2; exit 1; }
                prefix="$2"; shift ;;
            --dry-run)      dry=1 ;;
            --delete-calcs) delete=1 ;;
            --force)        force=1 ;;
            --gpu)
                [ $# -ge 2 ] || { echo "--gpu needs a GPU type, e.g. --gpu rtx5080" >&2; exit 1; }
                gpu="$2"; shift ;;
            --mem)
                [ $# -ge 2 ] || { echo "--mem needs a size, e.g. --mem 4G" >&2; exit 1; }
                check_mem "$2"
                mem="$2"; shift ;;
            -*) echo "Unknown option: $1" >&2; exit 1 ;;
            *)  inputs+=("$(local_path "$1")") ;;
        esac
        shift
    done
    # No files or folders given: the configs for run_all, relative to the repo root.
    [ ${#inputs[@]} -gt 0 ] || inputs=("$(pwd)/comp_configurations")

    # The config files: every *.json of a folder (sorted by name), or the file itself.
    files=()
    for p in "${inputs[@]}"; do
        if [ -d "$p" ]; then
            found=0
            for f in "${p%/}"/*.json; do
                [ -f "$f" ] || continue
                files+=("$f"); found=1
            done
            [ "$found" = 1 ] || { echo "There are no *.json files in $p." >&2; exit 1; }
        elif [ -f "$p" ]; then
            files+=("$p")
        else
            echo "$p does not exist (give config files or folders with them)." >&2
            exit 1
        fi
    done

    # Run name = prefix + file name without .json; valid, distinct, and the files are JSON.
    names=()
    declare -A source_of=()
    for f in "${files[@]}"; do
        base=$(basename "$f")
        name="$prefix${base%.json}"
        check_name "$name"
        if [ -n "${source_of[$name]:-}" ]; then
            echo "Two configs give the same run name '$name':" >&2
            echo "    ${source_of[$name]}" >&2
            echo "    $f" >&2
            exit 1
        fi
        source_of[$name]="$f"
        check_config_file "$f"
        names+=("$name")
    done

    set_gres "$gpu"

    # Names already used on the cluster (finished, running or queued), all in one ssh call.
    # The names contain only safe characters (check_name), so they can go into the command.
    taken=$(rssh "cd '$REMOTE_DIR' 2>/dev/null || exit 0; for n in ${names[*]}; do { [ -e \"results/\$n\" ] || [ -e \"jobs/\$n\" ]; } && echo \"\$n\"; done; true")
    if [ -n "$taken" ] && [ "$force" != 1 ]; then
        echo "These run names are already used on $CLUSTER_HOST (finished or queued):" >&2
        printf '    %s\n' $taken >&2
        echo "Pick other names (e.g. --prefix), or use --force to delete the old runs and start again." >&2
        exit 1
    fi

    commit=$(local_commit)
    echo ">> ${#names[@]} run(s) (commit $commit, GPU ${gpu:-any}${mem:+, host RAM $mem}${taken:+; --force deletes the old runs of the taken names first}):"
    for name in "${names[@]}"; do
        mark=""; printf '%s\n' "$taken" | grep -qxF -- "$name" && mark="   (taken: the old run is deleted)"
        printf '    %-30s <- %s%s\n' "$name" "${source_of[$name]}" "$mark"
    done
    if [ "$dry" = 1 ]; then
        echo ">> Dry run: nothing was deleted, uploaded or submitted."
        exit 0
    fi

    if [ -n "$taken" ]; then
        # --force does NOT stop a job that is still running under such a name (see the top).
        echo ">> Deleting the old results and uploads of the taken names (--force)"
        rssh "cd '$REMOTE_DIR' && for n in $(printf '%s ' $taken); do rm -rf \"results/\$n\" \"jobs/\$n\"; done"
    fi

    # The settings sbatch_run reads; the same as a plain submit.
    mode=submit; base_short=""; new_short=""; repeat=""; bench_set=""
    extra=(); [ -n "$mem" ] && extra+=(--mem="$mem")
    submitted=()
    for name in "${names[@]}"; do
        if upload_sources "$name" "${source_of[$name]}" && out=$(sbatch_run "$name"); then
            echo "$out"
            submitted+=("$name ${out##* }")
        else
            # A failed upload can leave jobs/<name> behind, which would keep the name taken
            # (sbatch_run deletes it itself when sbatch fails). If even this fails, `clean` does it.
            rssh "rm -rf '$REMOTE_DIR/jobs/$name'" || true
            echo "Submitting '$name' failed. Not submitted:" >&2
            started=0
            for n in "${names[@]}"; do
                [ "$n" = "$name" ] && started=1
                [ "$started" = 1 ] && echo "    $n" >&2
            done
            echo "Fix the problem and submit them again (the names are still free)." >&2
            break
        fi
    done
    echo ">> Submitted ${#submitted[@]} of ${#names[@]} run(s):"
    for s in ${submitted[@]+"${submitted[@]}"}; do
        printf '    %-30s job %s\n' "${s% *}" "${s##* }"
    done
    [ ${#submitted[@]} -eq ${#names[@]} ]
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
# Progress of the regression test or benchmark in results/$2 (job $1), from the progress file
# that tests/regression.sh and tests/benchmark.sh rewrite at every step. No estimate of the time
# left: the cases take very different times.
test_progress() {
    local phase a b c d e kind="regression test"
    [ "$(cat "results/$2/info/setup/mode.txt" 2>/dev/null)" = bench ] && kind=benchmark
    phase=""
    [ -f "results/$2/progress" ] && read -r phase a b c d e < "results/$2/progress"
    case "$phase" in
        build)   echo "job $1: $kind, compiling the $a version" ;;
        compare) echo "job $1: $kind, comparing the results" ;;
        run)
            # a = this run, b = number of runs, c = start of this run, d = case, e = version.
            echo "job $1: $kind, run $a/$b ($d, $e version), running for $(duration $(( $(date +%s) - c )))" ;;
        # Tests submitted before regression.sh wrote the progress file, or the job is still
        # unpacking the sources.
        *)       echo "job $1: $kind running (see slurm-$1.out)" ;;
    esac
}

squeue -u "$USER"
squeue -u "$USER" -h -t RUNNING -o '%i %j' | while read -r id name; do
    # A regression test or benchmark runs many short cases, so it reports which run it is at.
    if [ -f "results/$name/info/setup/baseline.tar.gz" ]; then
        test_progress "$id" "$name"
        continue
    fi
    if grep -q '^Computation time:' "results/$name/info/setup/parameters.txt" 2>/dev/null; then
        echo "job $id: computation finished, now plotting"
        continue
    fi
    failure=$(sed -n 's/^Run failed: *//p' "results/$name/info/setup/parameters.txt" 2>/dev/null)
    if [ -n "$failure" ]; then
        echo "job $id: computation FAILED ($failure), now plotting the frames"
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
    if [ -n "$2" ] && grep -q '^Run failed:' "$2" 2>/dev/null; then
        echo "FAILED ($(sed -n 's/^Run failed: *//p' "$2")), now plotting"; return
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

    date="-"; comp_time=""; failure=""; d=""
    if [ -n "$params" ]; then
        d=$(sed -n 's/^Run date: *//p' "$params" | cut -c1-16)
        [ -n "$d" ] && date=$d
        comp_time=$(sed -n 's/^Computation time: *//p' "$params" | sed 's/ *(.*//')
        failure=$(sed -n 's/^Run failed: *//p' "$params")
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
    # A regression test or benchmark (cluster.sh test or bench) keeps the baseline next to the
    # sources; it has no frames of its own (they are in new/ and base/), and report.txt is
    # written at its end. mode.txt tells the two apart (tests from before it are regression tests).
    is_test=0; kind="regression test"
    [ -f "$r/info/setup/baseline.tar.gz" ] && is_test=1
    [ "$(cat "$r/info/setup/mode.txt" 2>/dev/null)" = bench ] && kind=benchmark

    case "$jstate" in
        PENDING) state=queued;  details="job $jid" ;;
        RUNNING)
            state=running
            if [ "$is_test" = 1 ]; then details="job $jid on $node: $kind"
            else details="job $jid on $node: $(progress "$jid" "$params")"; fi ;;
        "")
            if [ "$is_test" = 1 ]; then
                if [ -f "$r/report.txt" ]; then
                    state=tested; details="$(tail -n 1 "$r/report.txt")"
                else
                    state=stopped
                    log=$(find_log "$name")
                    details="$kind did not finish${log:+; see $log}"
                fi
            elif [ ! -d "$r" ]; then
                # A job that failed before it moved its upload to the results left this behind.
                log=$(find_log "$name")
                state=leftover; details="upload in jobs/ without a job${log:+ (it failed at the start, see $log)}; run clean"
            elif [ -n "$failure" ]; then
                # The solver stopped itself: the solution broke down (the frames up to then are plotted)
                # or a frame could not be written.
                state=failed; details=$failure
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
    # Everything needed to decide how to submit: per GPU type the cards that are free and the ones
    # a job of ours could really get now (its node also has the CPUs and RAM free), the GPUs that
    # pending jobs of all users wait for, each free card with what blocks it, your own jobs, and
    # a verdict. "Our job" is what run.sbatch asks for, or --mem SIZE as with submit.
    # The remote script is sent on stdin; its arguments are the CPUs and MB of RAM we need.
    our_job_request "$@"
    "$SSH" "$CLUSTER_HOST" bash -s -- "$need_cpu" "$(mem_to_mb "$mem")" <<'REMOTE_SCRIPT'
{
    # "node NAME STATE CPUS(A/I/O/T) MEMORY_MB ALLOCMEM_MB gpu:TYPE:COUNT(S:0) gpu:TYPE:USED(IDX:..)"
    sinfo -N -h -p gp -O "NodeList:20,StateLong:20,CPUsState:30,Memory:20,AllocMem:20,Gres:80,GresUsed:80" | sed 's/^/node /'
    # Every pending array task on its own line (-r): "pending USER gres/gpu:TYPE:N", "gres/gpu:N" or "N/A"
    squeue -r -h -t PENDING -O "UserName:40,tres-per-node:80" | sed 's/^/pending /'
    # Our own jobs: "mine ID|NAME|STATE|REASON|START|ELAPSED|NODES"
    squeue -h -u "$USER" -o "%i|%j|%T|%r|%S|%M|%N" | sed 's/^/mine /'
} | awk -v need_cpu="$1" -v need_mem="$2" '
$1 == "node" && $7 ~ /^gpu:/ {
    n = $2; state = $3
    split($7, g, ":"); t = g[2]; total = g[3] + 0                # "gpu:rtx3060:2(S:0)"   -> 2
    split($8, u, ":"); free = total - (u[3] + 0)                 # "gpu:rtx3060:1(IDX:0)" -> 1 used
    split($4, c, "/"); cfree = c[2]; mfree = $5 - $6
    tot[t] += total; usd[t] += total - free
    nodes[t] = nodes[t] (nodes[t] == "" ? "" : " ") n
    if (state ~ /drain|down|fail|maint|reserv|not_resp|unk|inval/) {    # cannot get new jobs
        down[t] += free; nodes[t] = nodes[t] "(" state ")"; next
    }
    if (free < 1) next
    fre[t] += free
    why = ""
    if (cfree < need_cpu) why = why ", too few CPUs"
    if (mfree < need_mem) why = why ", too little RAM"
    line = sprintf("  %-5s %-8s %d card(s)   CPU %2d free   RAM %5.1f GB free   -> ", n, t, free, cfree, mfree / 1024)
    if (why == "") {
        line = line "FITS"; usable[t] += free; fits = fits " " n "(" t ")"
    } else {
        line = line "no" why
        if (cfree >= need_cpu && mfree >= 2048) {     # less would not be enough for our jobs
            line = line sprintf(" (fits with --mem %dG)", int(mfree / 1024)); smaller = 1
        }
    }
    if (state ~ /-$/) line = line ", planned for a waiting job"
    freelines[++nf] = line
}
$1 == "pending" && $3 ~ /gpu/ {
    k = split($3, p, ":"); cards = p[k] + 0; if (cards < 1) cards = 1
    if (k >= 3) wait[p[2]] += cards; else wait_any += cards
    waiting += cards
}
$1 == "mine" {
    split(substr($0, 6), m, "|")
    if (m[3] == "PENDING") s = sprintf("PENDING (%s), Slurm estimates the start at %s", m[4], m[5])
    else                   s = sprintf("%s on %s for %s", m[3], m[7], m[6])
    mine[++nm] = sprintf("  %-8s %-20s %s", m[1], m[2], s)
}
END {
    printf "%-9s %5s %7s %5s %5s %7s %8s   %s\n", "TYPE", "CARDS", "IN USE", "DOWN", "FREE", "USABLE", "WAITING", "NODES"
    cmd = "sort"
    for (t in tot) {
        printf "%-9s %5d %7d %5d %5d %7d %8d   %s\n", t, tot[t], usd[t], down[t], fre[t], usable[t], wait[t], nodes[t] | cmd
        all_usable += usable[t]
    }
    close(cmd)
    printf "%-9s %5s %7s %5s %5s %7s %8d\n", "any type", "", "", "", "", "", wait_any
    print "DOWN = idle cards on drained or down nodes. FREE = idle cards on working nodes."
    print "USABLE = free cards whose node also has the CPUs and RAM our job needs."
    print "WAITING = GPUs wanted by the pending jobs of all users (array tasks counted one by one)."
    print ""
    print "Free cards:"
    if (nf == 0) print "  none"
    for (i = 1; i <= nf; i++) print freelines[i]
    print ""
    print "Your jobs:"
    if (nm == 0) print "  none queued or running"
    for (i = 1; i <= nm; i++) print mine[i]
    print ""
    if (all_usable > 0) {
        print "-> A job submitted now can start at once on:" fits
        if (waiting > 0) printf "   (%d GPU(s) are wanted by waiting jobs; one with a higher priority may get it first.)\n", waiting
        print "   Any type: cluster.cmd submit NAME. Only one type: add --gpu TYPE."
    } else {
        print "-> No free card is usable now: a job submitted now has to wait."
        if (smaller) print "   A smaller RAM request would fit somewhere (see above); our jobs have peaked at 1.7 GB."
        print "   After submitting, cluster.cmd status shows why it waits and Slurm'"'"'s estimate of the start."
    }
    print "   Who holds the nodes: cluster.cmd nodes."
}'
REMOTE_SCRIPT
    ;;
nodes)
    # Every node of the gp partition: free GPUs, CPUs and RAM, the running jobs grouped by user
    # and job name, and whether a job of ours would fit there now. Free RAM is RAM that no job
    # has reserved (Slurm schedules by the reservations, not by what is really in use).
    # "Our job" is what run.sbatch asks for (1 GPU, its CPUs and RAM), or --mem SIZE as with submit.
    our_job_request "$@"
    "$SSH" "$CLUSTER_HOST" bash -s -- "$need_cpu" "$(mem_to_mb "$mem")" <<'REMOTE_SCRIPT'
{
    # "node NAME STATE CPUS(A/I/O/T) MEMORY_MB ALLOCMEM_MB GRES GRES_USED"
    sinfo -N -h -p gp -O "NodeList:20,StateLong:20,CPUsState:30,Memory:20,AllocMem:20,Gres:80,GresUsed:80" | sed 's/^/node /'
    # "job NODELIST USER cpu=2,mem=8G,node=1,billing=8,gres/gpu=1 NAME" (the name last: it may contain spaces)
    squeue -h -t RUNNING -O "NodeList:40,UserName:30,tres-alloc:200,Name:200" | sed 's/^/job /'
} | awk -v need_cpu="$1" -v need_mem="$2" '
function mb(s,   v) {                                     # "8G" -> 8192, "1500M" -> 1500
    v = s + 0
    if (s ~ /T$/) v *= 1048576; else if (s ~ /G$/) v *= 1024; else if (s ~ /K$/) v /= 1024
    return v
}
$1 == "node" {
    n = $2; order[++nn] = n; state[n] = $3
    split($4, c, "/"); cfree[n] = c[2]; ctot[n] = c[4]
    mtot[n] = $5; mfree[n] = $5 - $6
    gtot[n] = 0; gfree[n] = 0; gtype[n] = "none"
    if ($7 ~ /^gpu:/) {
        split($7, g, ":"); gtype[n] = g[2]; gtot[n] = g[3] + 0      # "gpu:rtx3060:2(S:0)"   -> 2
        split($8, u, ":"); gfree[n] = gtot[n] - (u[3] + 0)          # "gpu:rtx3060:1(IDX:0)" -> 1 used
    }
}
$1 == "job" {
    name = $5; for (i = 6; i <= NF; i++) name = name " " $i
    cpu = 0; mem = 0; gpu = 0
    k = split($4, t, ",")
    for (i = 1; i <= k; i++) {
        split(t[i], kv, "=")
        if (kv[1] == "cpu") cpu = kv[2]
        else if (kv[1] == "mem") mem = mb(kv[2])
        else if (kv[1] == "gres/gpu") gpu = kv[2]
    }
    key = $2 SUBSEP $3 SUBSEP name
    if (!(key in cnt)) { nk[$2]++; keys[$2, nk[$2]] = key }
    cnt[key]++; kcpu[key] += cpu; kmem[key] += mem; kgpu[key] += gpu
}
END {
    for (i = 1; i <= nn; i++) {
        n = order[i]
        why = ""
        if (state[n] ~ /drain|down|fail|maint|not_resp|unk/) why = why ", node " state[n]
        if (gfree[n] < 1)          why = why ", no free GPU"
        if (cfree[n] < need_cpu)   why = why ", too few free CPUs"
        if (mfree[n] < need_mem)   why = why ", too little free RAM"
        fits = (why == "") ? "FITS" : "no" why
        printf "%-5s %-10s GPU %d/%d free (%s)  CPU %d/%d free  RAM %.1f/%.0f GB free  -> %s\n", \
            n, state[n], gfree[n], gtot[n], gtype[n], cfree[n], ctot[n], mfree[n] / 1024, mtot[n] / 1024, fits
        for (j = 1; j <= nk[n]; j++) {
            key = keys[n, j]; split(key, f, SUBSEP)
            printf "      %-10s %3d job(s) %4d CPU %6.1f GB %2d GPU   %s\n", \
                f[2], cnt[key], kcpu[key], kmem[key] / 1024, kgpu[key], f[3]
        }
    }
    print ""
    print "Free RAM = RAM that no job has reserved; Slurm schedules by the reservations, not by real use."
    print "A state ending in \"-\" (e.g. mixed-) means the node is planned for a pending job with a higher"
    print "priority; a job of ours can still start there if its time limit ends before that job is due."
    print "Your pending jobs and why they wait: cluster.cmd status. Other RAM request: --mem SIZE."
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
    sed -n '2,/^$/p' "$0"                 # the comment at the top, up to the first empty line
    exit 1
    ;;
esac
