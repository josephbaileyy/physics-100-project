# ph100-delta-scutis

Lomb-Scargle period analysis for PH 100 variable-star targets.

The current primary workflow is the FITS reduction and differential
photometry analysis for V0399 UMa and AM CVn. See `docs/ANALYSIS_PLAN.md` for the
full reduction plan, calibration caveats, and acceptance checks. See
`docs/RESULTS.md` for the current V0399 and AM CVn results, uncertainty tests,
and model comparison summary. `docs/FINAL_PROJECT_REPORT_DRAFT.tex` and
`docs/FINAL_PROJECT_REPORT_DRAFT.pdf` contain a rough final report draft.

## Layout

- `docs/ANALYSIS_PLAN.md` — authoritative project workflow and QC plan
- `docs/RESULTS.md` — current V0399 and AM CVn results and uncertainty summary
- `docs/FINAL_PROJECT_REPORT_DRAFT.md/.tex/.pdf` — rough final report draft
- `configs/targets.json` — raw data paths, aperture settings, and manual target/comparison-star coordinates
- `scripts/project_pipeline.py` — FITS inventory, calibration, reduction, source proposals, photometry, and periodograms
- `scripts/v0399_contamination_tests.py` — V0399 blending, difference-image, model, and injection-recovery tests
- `scripts/amcvn_sequence_analysis.py` — AM CVn AAVSO sequence photometry and period analysis
- `analysis/project/` — generated project CSV summaries and plots that are useful for sharing; FITS files are ignored
- `data/` — historical ASAS-SN/AAVSO light curves
- `lomb_scargle.py` — standalone CSV Lomb-Scargle helper
- `results/` — historical standalone CSV periodogram outputs
- `PracticalLombScargle/` — VanderPlas (2018) reference paper as a submodule

## Project FITS Workflow

The repository is prepared to share code, reports, CSV summaries, and plots.
Raw FITS data and reduced FITS products are intentionally ignored by Git because
they are several GB. To reproduce the full reduction, obtain the raw data from
the shared PHYSICS100 Google Drive paths described in `docs/ANALYSIS_PLAN.md`.

```bash
pip install -r requirements.txt

python scripts/project_pipeline.py inventory
python scripts/project_pipeline.py calibrate
python scripts/project_pipeline.py reduce --target V0399_UMa
python scripts/project_pipeline.py reduce --target AM_CVn
python scripts/project_pipeline.py propose-stars --target V0399_UMa
python scripts/project_pipeline.py propose-stars --target AM_CVn
```

Inspect the proposal tables and overlays under `analysis/project/qc/`, then fill
the `x`/`y` values in `configs/targets.json` for the target, comparison,
and check stars. After that:

```bash
python scripts/project_pipeline.py photometry --target V0399_UMa
python scripts/project_pipeline.py photometry --target AM_CVn
python scripts/project_pipeline.py periodogram --target V0399_UMa
python scripts/project_pipeline.py periodogram --target AM_CVn
```

`run-all` performs the whole reproducible pass and stops after source proposals
if manual coordinates have not been filled in yet.

Additional final-analysis scripts:

```bash
python scripts/v0399_contamination_tests.py
python scripts/v0399_polish_figures.py
python scripts/amcvn_sequence_analysis.py
```

The most useful shareable products are:

- `docs/RESULTS.md`
- `docs/FINAL_PROJECT_REPORT_DRAFT.pdf`
- `analysis/project/contamination_tests/V0399_UMa/plots/V0399_final_summary_panel.png`
- `analysis/project/plots/AM_CVn_sequence_X42421BZ_light_curve.png`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/short_5_to_90min_periodogram.png`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/AM_CVn_sequence_reference_overlay.png`

## Historical CSV Usage

```bash
pip install -r requirements.txt
python lomb_scargle.py "data/ASASSN-V_J093023.98+612034.3.csv"
python lomb_scargle.py "data/ASASSN-V_J123454.75+373746.7.csv"
```

The script applies Lomb-Scargle to every valid observation on the original HJD
timeline. It does not condense the light curve into one or two phase periods
before fitting; the two-cycle folded plot is produced only afterward to make the
best period easier to inspect visually.

## Result so far

`ASASSN-V J093023.98+612034.3` → P ≈ 0.0615 d (88.6 min), FAP ≈ 4×10⁻⁷². Folded light curve shows an asymmetric ~0.5 mag pulsation consistent with a δ Scuti.

`ASASSN-V J123454.75+373746.7` → strongest unconstrained peak at P ≈
0.9984 d, FAP ≈ 6.7×10⁻¹⁰. Because that is very close to the one-day
ground-based observing cadence, also inspect shorter-period searches such as
`--max-period 0.3`; that gives P ≈ 0.1998 d, FAP ≈ 2.5×10⁻⁴.
