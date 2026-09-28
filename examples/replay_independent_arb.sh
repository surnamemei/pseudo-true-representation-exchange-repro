#!/usr/bin/env sh
# Independent Arb replays: z = 2 crossing certificates (R1); continuum partition for N >= 35,377 (R4, after R3).
set -e
python replay/replay_and_compare.py all R1 R3 R4
