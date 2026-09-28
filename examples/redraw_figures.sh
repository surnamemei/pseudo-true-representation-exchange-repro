#!/usr/bin/env sh
# Redraw the manuscript figures from the archived numerical results (no new search).
set -e
python replay/replay_and_compare.py figures
