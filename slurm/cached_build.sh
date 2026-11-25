#!/bin/bash
# Builds the solver in the source tree SRC_DIR, or takes the binary from the build cache when
# the same code was already built for this GPU type, and prints the path of the binary on
# stdout (everything else goes to stderr). Run on the compute node, with CUDA on the PATH.
#
#   bash slurm/cached_build.sh SRC_DIR CACHE_DIR [COMMIT]
#
# COMMIT only labels a new cache entry (info.txt). The build uses `make -j N -O` with N =
# SLURM_CPUS_PER_TASK (default 2): -O keeps the messages of each file together.
#
# Cache entry: CACHE_DIR/sm_XX-KEY/ with Phase-field, info.txt (where and from what it was built)
# and last_used (the time of the last use in seconds since 1970, which `cluster.sh cache prune`
# goes by; in the content, not the file's time: a `touch` from a compute node never reaches gp3
# over NFS, tested 2026-11-20, while a new file renamed into place does). KEY is a hash of
# everything the binary depends on:
#   - every file of SRC_DIR except the ones that cannot affect the build (config/, python/,
#     tests/, slurm/, results, build/, comp_configurations/, databases/, *.md, .git*, requirements.txt).
#     The rest counts even if the Makefile does not use it (old commits had their sources in the
#     top folder), so an unknown file can only cause a needless build, never a wrong binary;
#   - the GPU's compute capability (the Makefile compiles with -arch=native);
#   - nvcc --version, the gcc version, and the contents of the system TNL (/usr/include/TNL).
# Two jobs with the same code may both build and store it: the binary is copied under a
# temporary name and renamed, so a job sees no binary or a complete one, and the last one wins
# (both are built from the same inputs). A cache that cannot be read or written never fails
# the job: it builds, or keeps its own binary.

set -euo pipefail

src=$(cd "$1" && pwd)
cache=$2
commit=${3:-unknown}
jobs=${SLURM_CPUS_PER_TASK:-2}

log() { echo "Build cache: $*" >&2; }

# Builds in $src and sets $built to the binary. Plain `make` builds the binary in every version
# of the Makefile (it moved from build/ to build/sm_XX/ over time).
build() {
    echo ">> Building $src (make -j $jobs)" >&2
    make -C "$src" -j "$jobs" -O --no-print-directory >&2 || { echo "The build of $src failed" >&2; exit 1; }
    built=$(find "$src/build" -name Phase-field -type f -perm -u+x | sed -n 1p)
    [ -n "$built" ] || { echo "No Phase-field binary in $src/build after the build" >&2; exit 1; }
}

# The same command as GPU_TAG in the Makefile (sed instead of head: head can end the pipe early,
# which fails it under pipefail).
cc=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | sed -n 1p | tr -d '. \r') || cc=""
if [ -z "$cc" ]; then
    log "no GPU found, building without the cache"
    build
    echo "$built"
    exit 0
fi

code_files() {
    cd "$src"
    find . \( -path ./config -o -path ./python -o -path ./tests -o -path ./slurm -o -path ./results \
              -o -path ./build -o -path ./comp_configurations -o -path ./databases \) -prune \
         -o -type f ! -name '*.md' ! -name '.git*' ! -name requirements.txt -print0 \
        | LC_ALL=C sort -z | xargs -0 -r sha256sum
}
key=$( {
    echo "cached_build.sh key, version 1"
    echo "sm_$cc"
    nvcc --version
    gcc --version | sed -n 1p
    find /usr/include/TNL -type f -print0 | LC_ALL=C sort -z | xargs -0 -r sha256sum
    code_files
} | sha256sum | cut -c 1-16)
name="sm_$cc-$key"
entry="$cache/$name"
# Temporary files of this job in the entry (renamed into place when complete).
tmp="$entry/.tmp-${SLURM_JOB_ID:-$$}-${HOSTNAME%%.*}"

mark_used() {
    date +%s > "$tmp.used" && mv -f "$tmp.used" "$entry/last_used"
}

# Hit: copy the binary to the node's disk, so pruning the cache cannot touch a running job.
local_bin="$src/build/cached/Phase-field"
if [ -f "$entry/Phase-field" ]; then
    if mkdir -p "$(dirname "$local_bin")" && cp "$entry/Phase-field" "$local_bin" && chmod u+x "$local_bin"; then
        mark_used 2>/dev/null || rm -f "$tmp.used" 2>/dev/null || true
        log "using $name, no build ($(sed -n 1p "$entry/info.txt" 2>/dev/null || echo 'no info.txt'))"
        echo "$local_bin"
        exit 0
    fi
    log "could not copy $entry/Phase-field, building instead"
fi

log "nothing for $name yet, building"
build

# Store. info.txt and last_used come first, so an entry with a binary always has them.
store() {
    mkdir -p "$entry" &&
    {
        echo "Built $(date '+%Y-%m-%d %H:%M') on ${HOSTNAME%%.*} ($(nvidia-smi --query-gpu=name --format=csv,noheader | sed -n 1p)) by job ${SLURM_JOB_ID:-none}${FOLDER:+ (run $FOLDER)}, commit $commit"
        echo "nvcc: $(nvcc --version | sed -n 's/.*release /release /p')"
        echo "gcc: $(gcc --version | sed -n 1p)"
    } > "$tmp.info" &&
    cp "$built" "$tmp.bin" &&
    mv -f "$tmp.info" "$entry/info.txt" &&
    mark_used &&
    mv -f "$tmp.bin" "$entry/Phase-field"
}
if store; then
    log "stored as $name"
else
    rm -f "$tmp.info" "$tmp.bin" "$tmp.used" 2>/dev/null || true
    log "could not store the binary in $entry (the run goes on)"
fi
echo "$built"
