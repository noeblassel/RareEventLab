#!/usr/bin/env bash
# Install Julia, the notebook environment, and a Jupyter kernel, all inside this folder:
#   julia/         the Julia binaries
#   .julia-depot/  the Julia packages (and Jupyter, if IJulia has to install it)
# Nothing is written outside this folder, except the small Jupyter kernel file.
# Uninstall with ./uninstall.sh, or delete this folder and run
#   jupyter kernelspec remove julia-rareeventlab
# Usage: ./setup.sh [number_of_threads]      (default: number of cores)
# Linux and macOS. On Windows, run it inside WSL.
set -euo pipefail

cd "$(dirname "$0")"
HERE="$(pwd)"
THREADS="${1:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu)}"
VERSION="1.13.3"
MINOR="${VERSION%.*}"

case "$(uname -s)-$(uname -m)" in
    Linux-x86_64)               URL="linux/x64/$MINOR/julia-$VERSION-linux-x86_64.tar.gz" ;;
    Linux-aarch64)              URL="linux/aarch64/$MINOR/julia-$VERSION-linux-aarch64.tar.gz" ;;
    Darwin-arm64)               URL="mac/aarch64/$MINOR/julia-$VERSION-macaarch64.tar.gz" ;;
    Darwin-x86_64)              URL="mac/x64/$MINOR/julia-$VERSION-mac64.tar.gz" ;;
    *) echo "Unsupported platform: $(uname -s) $(uname -m)" >&2; exit 1 ;;
esac

if [ ! -x julia/bin/julia ]; then
    echo ">> Downloading Julia $VERSION into $HERE/julia"
    rm -rf julia
    mkdir julia
    curl -fL --progress-bar "https://julialang-s3.julialang.org/bin/$URL" |
        tar -xz -C julia --strip-components 1
fi

JULIA="$HERE/julia/bin/julia"
# the packages go in .julia-depot; julia/share/julia holds the precompiled standard library
export JULIA_DEPOT_PATH="$HERE/.julia-depot:$HERE/julia/share/julia"
echo ">> Julia: $("$JULIA" --version)"

echo ">> Installing the notebook environment (this takes a few minutes)"
"$JULIA" --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo ">> Installing the Jupyter kernel 'Julia RareEventLab' with $THREADS threads"
"$JULIA" --project=. -e "using IJulia; IJulia.installkernel(\"Julia RareEventLab\", \"--project=$HERE\";
    specname = \"julia-rareeventlab\",
    env = Dict(\"JULIA_NUM_THREADS\" => \"$THREADS\", \"JULIA_DEPOT_PATH\" => \"$JULIA_DEPOT_PATH\"))"

if ! command -v jupyter >/dev/null 2>&1; then
    echo ">> No 'jupyter' found: you can start one from Julia with"
    echo "     JULIA_DEPOT_PATH=\"$JULIA_DEPOT_PATH\" ./julia/bin/julia --project=. -e 'using IJulia; notebook(dir = pwd())'"
fi
echo ">> Done. Open notebook.ipynb with the kernel 'Julia RareEventLab'."
