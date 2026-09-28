# Code and replay guide for the anonymous paper

## Environments

- Main paper and supplement: an existing WSL TeX Live/latexmk environment.
- Original interval certificates: Python with mpmath and the dependencies already listed in REPRODUCIBILITY_README.md.
- Independent all-large-N replay: Linux or WSL Python 3.10+ with python-flint 0.9.0. It imports no Stage 14/15 checker functions.
- Independent asymptotic audit: Python with SymPy and mpmath.

## Shortest Stage 16 replay

From the project root:

    python stage16_symbolic_asymptotics.py

In WSL after installing python-flint into an isolated environment:

    python3 -m venv /tmp/finite-tone-stage16-venv
    /tmp/finite-tone-stage16-venv/bin/pip install python-flint==0.9.0
    /tmp/finite-tone-stage16-venv/bin/python stage16_independent_arb.py

The Arb script prints progress and writes independent_N0_replay.json and independent_N0_replay.csv. Expect 27,708 passed frozen cells, ten original cells refined (26 subdivision nodes), zero failures, and verified=true. The JSON includes directed margins. Running the scripts overwrites only their Stage 16 output files.

## Derived coefficient tables and PDFs

    python stage14_asymptotics.py
    python stage15_asymptotic_data.py

The first script performs analytical differentiation and numerical evaluation; it uses multiprecision before every division. The second produces the nine-length validation table and largeN_theory_vs_data.pdf. The direct N11/15/21/31/41 lambda values are read from existing certificates; N101/201/501/1001 are numerical solves. Neither data table is an interval proof.

## Compile anonymous files

    cd paper && latexmk -pdf -interaction=nonstopmode -halt-on-error main_tsp_reviewfriendly.tex
    cd ../supplement && latexmk -pdf -interaction=nonstopmode -halt-on-error supplement.tex

For the complete certificate dependency graph and historical correction notes, use REPRODUCIBILITY_README.md. The archive index identifies the primary proof records; artifact_manifest.csv supplies frozen hashes.
