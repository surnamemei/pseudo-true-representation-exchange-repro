#!/usr/bin/env sh
# All seven certificate chains in parallel, then the comparison with the archived records.
set -e
python replay/replay_and_compare.py all
