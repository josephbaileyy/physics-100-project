# V0399 UMa and AM CVn Analysis Plan

This repository now treats the local CCD observations as the primary workflow for
V0399 UMa and AM CVn. The science product is differential flux light curves and
Lomb-Scargle periodograms, not absolute calibrated V magnitudes.

Current V0399 results and uncertainty-reduction tests are documented in
`docs/RESULTS.md`. This file remains the workflow and QC plan.

## Data Inventory

Raw science lights:

- V0399 UMa: 109 V-band `8s` frames in the science light directory.
- AM CVn: 140 V-band frames split across `35s`, `40s`, and `45s`.
- Extra `light_V_8s` frames are present and should be treated as separate
  context data unless they are intentionally assigned to a target later.

May 20 calibration frames:

- Bias: 10 frames. The first five were taken around `2026-05-21T03:24` UTC
  and may have dome-light contamination. The last five were taken around
  `2026-05-21T07:27` UTC and are the default science bias set.
- Darks: matched dark sets exist for `8s`, `35s`, `40s`, and `45s`, plus flat
  exposure darks for `0.3s`, `4s`, and `6s`.
- Flats: V flats exist at `0.3s`, `4s`, and `6s`. The `0.3s` flats are the
  known dome-light-problem set.

The active configuration is `observations/targets.json`. It points at local
ignored copies under `observations/raw/may20/`, stores the output path, aperture settings,
saturation threshold, target light globs, and the manual object positions used
for photometry.

## Workflow

Run the pipeline from the repository root:

```bash
python scripts/project_pipeline.py inventory
python scripts/project_pipeline.py calibrate
python scripts/project_pipeline.py reduce --target V0399_UMa
python scripts/project_pipeline.py reduce --target AM_CVn
python scripts/project_pipeline.py propose-stars --target V0399_UMa
python scripts/project_pipeline.py propose-stars --target AM_CVn
```

After inspecting the proposal tables and overlays, fill the `x` and `y` values
for the target, comparison, and check stars in `observations/targets.json`.
Then run:

```bash
python scripts/project_pipeline.py photometry --target V0399_UMa
python scripts/project_pipeline.py photometry --target AM_CVn
python scripts/project_pipeline.py periodogram --target V0399_UMa
python scripts/project_pipeline.py periodogram --target AM_CVn
```

For a full rerun after coordinates are filled:

```bash
python scripts/project_pipeline.py run-all
```

Outputs are written under `analysis/project/`:

- `inventory/`: frame counts and header summary.
- `calibration/`: master biases, darks, flat candidates, and flat QC table.
- `reduced/<target>/`: calibrated science frames for the selected flat.
- `qc/<target>/`: reference-frame overlays and proposed source tables.
- `photometry/<target>/`: aperture photometry, differential light curves,
  saturation flags, and RMS diagnostics.
- `periodograms/<target>/`: time series, Lomb-Scargle periodograms, folded
  light curves, and result tables.

## Calibration Strategy

The default science calibration uses the late master bias, matched exposure
darks, and the best flat candidate selected by QC. The pipeline builds these
flat candidates:

- `flat_0p3_late_bias`: `0.3s` V flats minus late master bias and `0.3s` dark.
- `flat_0p3_early_bias`: `0.3s` V flats minus early master bias and `0.3s`
  dark, to test whether shared dome-light contamination cancels better.
- `flat_4s_masked`: star-masked/rejected `4s` V flats using late bias and
  `4s` dark.
- `flat_6s_masked`: star-masked/rejected `6s` V flats using late bias and
  `6s` dark.
- `no_flat`: diagnostic only, never preferred unless all flats visibly fail.

The selection metric favors a smooth normalized flat with low robust scatter
and few high-contrast residuals. It is still a QC choice: inspect the generated
flat images and comparison/check-star scatter before trusting the final run.

## Photometry Strategy

All source fluxes are divided by exposure time before differential photometry,
which lets AM CVn combine `35s`, `40s`, and `45s` frames. The ensemble
comparison flux is the median of unsaturated, positive comparison-star flux
rates in each frame.

AM CVn should use the AAVSO comparison sequence, including stars 96 and 99.
If 96 or 99 saturates in one or two frames, only that star is excluded for those
frames. A frame is rejected only when the target is saturated, registration
fails, target flux is invalid, or fewer than the configured minimum comparison
stars remain.

V0399 UMa has no sequence in this plan. Use `propose-stars` to identify
unsaturated, isolated field stars, then manually enter target/comparison/check
pixel coordinates in the config. If the brightest star saturates, flag and
exclude that object only; keep the frame if the target and enough comparison
stars remain valid.

## Period Search

Lomb-Scargle uses the continuous FITS `JD` timeline. The default period ranges
are target-specific in the config:

- V0399 UMa: 5 to 180 minutes.
- AM CVn: 5 to 90 minutes.

Each target also gets a wider sensitivity run by expanding the upper period
bound to at least 360 minutes. The reported table includes the best period,
frequency, power, false-alarm probability, and number of retained points.

## Acceptance Checks

- Inventory frame counts match the expected May 20 science and calibration
  counts.
- Calibrated reference frames look sane for the selected flat and the no-flat
  diagnostic does not produce a contradictory conclusion.
- Target and check-star plots show that the comparison ensemble is not driving
  the signal.
- Saturation and registration flags are saved and explain any rejected point.
- Periodogram peaks remain stable when marginal saturated-comparison frames are
  removed.
