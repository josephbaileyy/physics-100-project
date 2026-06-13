# AM CVn Time-Series Photometry

This repository contains the PHYSICS 100 final project workflow for V-band
time-series photometry of AM CVn. The project reduces a 2026 May 20 observing
sequence from the Stanford Student Observatory, performs differential aperture
photometry with AAVSO comparison stars, and models the resulting light curve to
recover AM CVn's positive superhump and first harmonic.

## Final Products

- [Final report](products/final_report/AM_CVn_final_report.pdf)
- [Final presentation](products/final_presentation/AM_CVn_final_presentation.pdf)

The submitted analysis finds a single-sinusoid period of `524.35 +/- 3.3 s`,
matching the first harmonic of AM CVn's double-humped waveform. A free
double-wave model gives the strongest independent period estimate,
`1050.05 +/- 6.7 s`, consistent with the published `1051.2 s` positive
superhump. The 111.12 minute baseline recovers the known signal, but it is too
short to improve the published period constraints or separate the nearby
orbital, positive-superhump, and negative-superhump clocks.

Older two-target and standalone CSV period-search work has been moved into
`archive/` so the active tree reflects the AM CVn final project.

## Repository Layout

- `src/amcvn/` - reusable analysis package for reduction helpers, photometry,
  uncertainty tests, model comparisons, and figure generation
- `scripts/` - small command-line entrypoints for the active AM CVn workflow
- `observations/targets.json` - active target configuration and local raw-data
  paths
- `analysis/am_cvn/` - generated AM CVn photometry tables, periodograms, model
  summaries, and final figures
- `analysis/calibration/` - calibration summaries and flat-field QC products
- `products/final_report/` - submitted final report PDF
- `products/final_presentation/` - submitted final presentation PDF
- `references/literature/am_cvn/` - bibliography and literature index
- `archive/` - retired exploratory work retained for provenance, but not part
  of the active public workflow

Raw and reduced FITS files are intentionally ignored because they are bulky and
environment-specific. The shareable outputs are the CSV summaries, PNG figures,
code, bibliography, report, and presentation.

## Reproducing The Workflow

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the FITS reduction pipeline:

```bash
python scripts/amcvn_pipeline.py inventory
python scripts/amcvn_pipeline.py calibrate
python scripts/amcvn_pipeline.py reduce --target AM_CVn
python scripts/amcvn_pipeline.py propose-stars --target AM_CVn
```

Inspect AM CVn proposal products under `analysis/am_cvn/qc/`, then fill any
needed manual positions in `observations/targets.json`. The sequence-specific
photometry and model analysis uses:

```bash
PYTHONPATH=src python -m amcvn.sequence_analysis
PYTHONPATH=src python -m amcvn.period_uncertainty
PYTHONPATH=src python -m amcvn.exploratory_period_modeling
PYTHONPATH=src python -m amcvn.physical_model_suite
PYTHONPATH=src python -m amcvn.stacked_model_fits
```

## Useful Outputs

- [AM CVn light curve](analysis/am_cvn/plots/AM_CVn_sequence_X42421BZ_light_curve.png)
- [Stacked model fits](analysis/am_cvn/plots/AM_CVn_stacked_model_fits.png)
- [Lomb-Scargle periodogram](analysis/am_cvn/periodograms/sequence_X42421BZ/short_5_to_90min_periodogram.png)
- [AAVSO comparison-star overlay](analysis/am_cvn/photometry/sequence_X42421BZ/AM_CVn_sequence_reference_overlay.png)
