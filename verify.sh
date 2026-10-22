#!/usr/bin/env bash
# Reproduce and check every numerical/symbolic result of the paper.   Usage: bash verify.sh [--quick|--full]
#   --quick (default): all generators, ~3 min CPU     --full: + finer time-domain grid (Nx=24000) and halved RK4 steps, ~25 min CPU
# Needs Python >= 3.10 with numpy>=2, scipy>=1.11, sympy, mpmath, matplotlib (pip install -r requirements.txt); no network, no compiler.
set -u
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
MODE="${1:---quick}"
case "$MODE" in --quick|--full) ;; *) echo "usage: bash verify.sh [--quick|--full]"; exit 2;; esac
"$PY" -c 'import sys; assert sys.version_info >= (3, 10), "Python >= 3.10 required"' || exit 2
"$PY" -c 'import numpy' 2>/dev/null || { echo "numpy is required: pip install -r requirements.txt"; exit 2; }
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
exec "$PY" verify.py "$MODE"
