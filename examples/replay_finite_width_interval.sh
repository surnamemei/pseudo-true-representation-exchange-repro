#!/usr/bin/env sh
# N = 21 finite-width amplitude interval (R6). The first command is stage17_epsilon_refine.py 0.035 14287.
set -e
python replay/replay_and_compare.py all R6
