# Sharing Notes

This repository is prepared for sharing code, reports, CSV summaries, and plots
with teammates.

## Included in Git

- Analysis scripts in `scripts/`
- Configuration in `configs/targets.json`
- Project documentation:
  - `docs/ANALYSIS_PLAN.md`
  - `docs/RESULTS.md`
  - `docs/FINAL_PROJECT_REPORT_DRAFT.md`
  - `docs/FINAL_PROJECT_REPORT_DRAFT.tex`
  - `docs/FINAL_PROJECT_REPORT_DRAFT.pdf`
- Shareable project products under `analysis/project/`, including:
  - CSV summaries and light curves
  - periodogram CSVs
  - final plots and QC overlays
  - AAVSO chart metadata for AM CVn sequence `X42421BZ`

## Excluded from Git

The raw and reduced FITS files are intentionally ignored:

- `raw/`
- `*.fit`, `*.fits`, `*.fts`

These files are several GB locally and should be shared through Google Drive or
another data archive, not ordinary GitHub. The scripts assume those raw data
paths are available locally if someone wants to reproduce the full reduction.

## Most Useful Files for Teammates

- `docs/RESULTS.md`
- `docs/FINAL_PROJECT_REPORT_DRAFT.pdf`
- `analysis/project/contamination_tests/V0399_UMa/plots/V0399_final_summary_panel.png`
- `analysis/project/plots/AM_CVn_sequence_X42421BZ_light_curve.png`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/short_5_to_90min_periodogram.png`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/AM_CVn_sequence_reference_overlay.png`
