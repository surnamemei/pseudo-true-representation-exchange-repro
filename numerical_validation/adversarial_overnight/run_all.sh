#!/usr/bin/env bash
# Launch the resumable adversarial validation runner, unbuffered, teeing to overnight.log.
# Safe inside tmux/screen; rerunning resumes from checkpoints.
set -o pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_CBWR=COMPATIBLE MPLBACKEND=Agg
export ADV_WORKERS="${ADV_WORKERS:-22}"
python -u run_all.py "$@" 2>&1 | tee -a overnight.log
