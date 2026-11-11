#!/usr/bin/env bash
# Remove what ./setup.sh installed: Julia, its packages, and the Jupyter kernel.
set -euo pipefail

cd "$(dirname "$0")"
rm -rf julia .julia-depot
if command -v jupyter >/dev/null 2>&1; then
    jupyter kernelspec remove -y julia-rareeventlab 2>/dev/null || true
fi
for d in "$HOME/.local/share/jupyter/kernels" "$HOME/Library/Jupyter/kernels"; do
    rm -rf "$d/julia-rareeventlab"
done
echo ">> Removed Julia, its packages and the 'Julia RareEventLab' kernel."
