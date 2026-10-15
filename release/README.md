# Release notes (draft; not published)

This directory prepares a later public snapshot. It does not create a remote and it does not push.

## What a reviewer would run

Script order, pinned packages, and the thread-cap proposal are in the files next to this note.

1. `release/script_order.md`
2. `release/requirements_freeze.txt`
3. `release/absolute_paths.md`
4. `release/cpu_cap.md`
5. `release/licenses.md`
6. `release/candidate_files.md`

The code is under the MIT License. Derived result tables and model files are under Creative Commons Attribution 4.0 International. Raw data matrices are not redistributed.

## What is not in git

`data/raw/`, `data/processed/` parquet files that are gitignored, and `.venv/` stay local. The non-imputed APOLLO text file is raw data and is not committed. `results/stage6/fake_null/` and `results/stage6/fake_planted/` were local checks and are not part of the release list.

## Figures

Final figures 3, 5, and 6 are `results/figures/fig3.pdf`, `fig5.pdf`, and `fig6.pdf`, with PNG at 300 dpi. English legends are `results/figures/legends_en.md`.
