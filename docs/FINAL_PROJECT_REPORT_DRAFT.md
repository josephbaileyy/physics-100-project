# Time-Series Photometry of V0399 UMa and AM CVn

**Rough draft for PHYSICS 100 final independent project report**  
**Team:** TODO  
**Authors:** TODO  
**Date:** TODO

## Abstract

We obtained time-series V-band observations of two suspected short-period
variable targets, V0399 UMa and AM CVn, on 2026 May 20 using the 0.7 m
telescope. Our goal was to reduce the CCD images into differential light curves
and test for periodic variability with Lomb-Scargle periodograms and fixed-period
model fits. Calibration was complicated by dome-light contamination in the
short flats and early bias frames, star contamination in the longer flats, and
tracking artifacts in some science images. We therefore compared multiple flat
fields, used per-frame source matching, and used control stars to quantify
systematic uncertainty. For V0399 UMa, the most consistent light curve has a
scatter near 1.4 percent, but no reliable detection of the cataloged 15.3 min
period. Difference-image and aperture-size tests indicate strong contamination
or common-mode coupling with nearby DN UMa. For AM CVn, AAVSO sequence
photometry gives a stronger signal: the dominant Lomb-Scargle peak is near
8.76 min with false-alarm probability 0.014, while fixed fits near 17.1 min are
weaker but still modestly favored. We conclude that AM CVn is detected as
variable, while the V0399 UMa data support only a cautious non-detection.

## 1. Introduction

Short-period stellar variability provides a compact way to probe stellar
structure, binary evolution, and accretion physics. In this project we observed
two targets selected for rapid variability: V0399 UMa and AM CVn. The scientific
goal was not to measure absolute calibrated magnitudes, but to test whether our
data could recover short-period signals in differential photometry.

AM CVn is the prototype of the AM CVn class: hydrogen-deficient ultracompact
interacting binaries in which a white dwarf accretes helium-rich material from a
compact donor. AM CVn systems have orbital periods of only minutes to about an
hour and are important laboratories for accretion physics and compact binary
evolution. They are also expected low-frequency gravitational-wave sources
(Levitan et al. 2013; Solheim 2010). The known AM CVn orbital timescale is about
17 min, so it is a natural target for a single-night high-cadence photometry
experiment.

V0399 UMa was included as a suspected rapid variable. AAVSO VSX lists a period
near 15.3 min for V0399 UMa. The field also contains DN UMa close to the target.
DN UMa is itself a variable object, so an important part of the analysis was to
test whether apparent V0399 UMa variability could be produced by blending,
tracking drift, or point-spread-function contamination from DN UMa.

Our specific project goals were:

1. Reduce the May 20 CCD images into calibrated science frames.
2. Produce differential flux light curves for V0399 UMa and AM CVn.
3. Search for short-period variability using Lomb-Scargle periodograms and
   fixed-period fits.
4. Quantify calibration, comparison-star, tracking, and blending systematics.
5. Decide which variability claims are supported by the data.

## 2. Observations

The observations were taken on 2026 May 20 in the V band. The V0399 UMa data
include two related image sets: 83 frames labeled `light_V_8s`, spanning
22.63 min, and 109 frames labeled `light_V0399_V_8s`, spanning 29.80 min. Both
sets were 8 s exposures. Together they span 54.11 min with a gap between the two
sequences. The `light_V_8s` frames were verified to contain V0399 UMa.

The AM CVn data include 140 V-band frames with exposure times of 35, 40, and
45 s. These frames span about 113 min. The exposure time changes were handled by
converting all aperture sums to flux rates before comparison. The final two
AM CVn science frames showed severe artifacts and very low peak counts, and
were excluded from the science light curve.

Calibration data included bias frames, dark frames matching the science
exposure times, and V-band flats with exposure times of 0.3, 4, and 6 s. The
calibration set had several known problems. All 0.3 s flats had dome-light
contamination. The 4 and 6 s flats contained stars. Bias frames taken early in
the night were affected by the dome-light issue, while later bias frames were
cleaner.

## 3. Data Reduction

Reduction was performed with Python using `astropy`, `sep`, `numpy`, `scipy`,
and `matplotlib`. We built a reproducible local pipeline that performs
calibration, source detection, aperture photometry, and period searches.

### 3.1 Calibration Frames

The default bias frame was the late master bias, because the early bias frames
were affected by dome-light contamination. Dark frames were median-combined by
exposure time and subtracted from science frames after bias subtraction.

Several flat-field candidates were compared:

- `flat_0p3_late_bias`: 0.3 s flats corrected with the late bias.
- `flat_0p3_early_bias`: 0.3 s flats corrected with the early contaminated bias.
- `flat_4s_masked`: star-masked 4 s flats, median combined.
- `flat_6s_masked`: star-masked 6 s flats, median combined.
- `flat_4s_masked_clipmean`: star-masked 4 s flats, sigma-clipped mean.
- `flat_6s_masked_clipmean`: star-masked 6 s flats, sigma-clipped mean.
- `no_flat`: a diagnostic reduction with no flat-field correction.

For the star-contaminated 4 and 6 s flats, we normalized each flat by its median,
subtracted a median-filtered smooth illumination model, masked compact positive
residuals, dilated the masks to include star wings, and combined the remaining
pixels. The best overall calibration choice for the science products was
`flat_4s_masked_clipmean`. This flat gave the lowest V0399 target scatter among
the tested flat fields, although the improvement over `no_flat` and the median
4 s flat was small.

Calibration uncertainty was quantified by comparing frame-by-frame normalized
V0399 differential flux across the real flat candidates. The combined median
calibration-induced scatter was 0.119 percent, with a 90th percentile of
0.353 percent. Including the `no_flat` diagnostic increased the median to
0.170 percent. Therefore calibration uncertainty was real, but smaller than the
final V0399 noise floor of about 1.4 percent.

### 3.2 Source Registration and Photometry

Simple fixed-aperture photometry did not work for V0399 UMa because tracking
drift shifted the sources significantly over time. We therefore used per-frame
source matching. For V0399 UMa, we detected V0399 and nearby DN UMa in each
frame, used the pair to define the local field transformation, and then matched
comparison stars to SEP detections in that frame. This fixed an early failure in
which comparison apertures drifted onto blank sky.

For AM CVn, the science images had no reliable WCS in the FITS headers, so we
used the AAVSO comparison sequence chart `X42421BZ`. The downloaded chart
contained labels 90, 96, 112, 125, 143, and 144. A good mid-run reference frame
was selected manually because the automatic minimum-ADU reference choice picked
one of the bad late artifact frames. The AAVSO sequence geometry was then fit
to detected source positions in the reference image.

Aperture photometry used circular apertures and sky annuli. For V0399 UMa, the
baseline aperture was 5 px with annulus 10--18 px. For AM CVn, an aperture grid
showed that 6 px with annulus 12--20 px gave the best target scatter among the
tested circular apertures.

For differential photometry, all fluxes were divided by exposure time to produce
flux rates. V0399 UMa used an ensemble of eight Gaia comparison stars common to
both V0399 pointings. AM CVn used the AAVSO sequence stars, but each comparison
star was first normalized by its own median flux before forming the comparison
ensemble. This step was necessary because the AAVSO sequence stars span
approximately V = 8.97 to V = 14.36; otherwise, the ensemble median can jump
when a bright saturated comparison star is excluded.

## 4. Analysis Methods

### 4.1 Lomb-Scargle Periodograms

We used Lomb-Scargle periodograms to search for sinusoidal periodic signals in
unevenly sampled time-series photometry. In the standard normalization, the
Lomb-Scargle power measures the fractional improvement of a sinusoidal model
over a constant model at each tested frequency. False-alarm probabilities were
computed using the `astropy.timeseries.LombScargle` implementation. These FAPs
must be interpreted carefully: they estimate the probability of obtaining a peak
at least this high under a noise-only null hypothesis, not the probability that
the star is non-variable.

For V0399 UMa, we searched periods from 3 to 45 min and also tested the catalog
period of 15.3 min directly. For AM CVn, we searched 5 to 90 min and a wider
3 to 180 min range. We also fit fixed sinusoids at 17.14 and 17.20 min because
those periods are close to the expected AM CVn orbital timescale.

### 4.2 Model Comparison and Controls

Because our light curves were limited by systematics, we did not rely only on
the highest Lomb-Scargle peak. We also compared constant models, fixed-period
sine models, sawtooth and eclipse-like templates, and two-star contamination
models using residual scatter and the Bayesian information criterion (BIC).

Control stars were essential. For V0399 UMa, DN UMa and HD 103187 were treated
as checks, not as guaranteed constants. For AM CVn, each AAVSO sequence star was
temporarily treated as a target and compared against the remaining sequence
stars to estimate the control-star scatter.

### 4.3 V0399 UMa Contamination Tests

DN UMa is near V0399 UMa and is itself variable. We therefore ran several tests
for blending and contamination:

- Aperture/annulus grid photometry over radii from 3 to 10 px.
- Difference-image residual photometry at V0399 UMa and DN UMa.
- Two-star mixture models in which V0399 and DN light curves could leak into
  each other.
- Injection-recovery tests for synthetic 15.3 min sine, dip, and sawtooth
  signals.
- Fixed-period tests at 15.3 min and nearby aliases.

These tests were designed to distinguish an intrinsic V0399 signal from
common-mode residuals, tracking effects, and contamination from DN UMa.

## 5. Results

### 5.1 Calibration and Data Quality

The calibration tests showed that no flat-field choice changed the qualitative
result for V0399 UMa. The best target-scatter calibration was
`flat_4s_masked_clipmean`, with V0399 scatter of 1.385 percent and DN/check
scatter of 1.673 percent. The no-flat diagnostic was surprisingly competitive,
with V0399 scatter of 1.394 percent. The dome-light early-bias correction for
the 0.3 s flats did not improve the science light curve.

For AM CVn, the most important data-quality issue was not the flat field, but
bad tracking/artifacts and the reference-frame choice. A contact sheet of early
and late frames showed that the final two AM CVn frames were unusable. These
frames were flagged as artifacts and removed. Several frames also saturated
comparison star 90; these points were excluded only for that star rather than
discarding the entire exposure.

![Figure 1. V0399 UMa robustness summary. The panels show the preferred
common-comparison light curve, aperture-size dependence, two-star/common-mode
model comparison, and injection recovery at the VSX period. The main result is
that V0399 UMa does not show a robust independent periodic signal; larger
apertures reduce scatter but increase V0399--DN correlation.](/Users/josephbailey/Documents/School/100/ph100-delta-scutis/analysis/project/contamination_tests/V0399_UMa/plots/V0399_final_summary_panel.png)

### 5.2 V0399 UMa

The most defensible V0399 UMa light curve uses the same eight Gaia comparison
stars in both V0399-related sequences. With this setup, the normalized
differential flux scatter is about 1.4 percent. The two V0399 sequences are
internally consistent when the common comparison ensemble is used; smaller
optimized comparison subsets sometimes reduce scatter, but they introduce
segment-dependent trends and are less reliable for combined physical claims.

The Lomb-Scargle periodogram did not recover a secure period. The best short
period in one optimized/balanced analysis was 4.68 min with amplitude about
0.40 percent, but the false-alarm probability was approximately 1.0. The
cataloged 15.3 min period produced only a weak fixed-period fit, with amplitude
about 0.28 percent and no convincing improvement over a constant model.

The contamination tests were decisive. As the aperture radius increased from
3 to 10 px, V0399 scatter decreased from about 3.65 percent to about
0.68 percent, but the V0399--DN correlation increased from about 0.17 to
0.78. Difference-image residuals at V0399 and DN were also highly correlated
with r about 0.916 for apertures of 4--6 px. This is not the signature of a
clean residual at V0399 alone. It is more consistent with shared subtraction
residuals, common-mode tracking effects, or blending between the nearby pair.

Injection-recovery tests show that a coherent 15.3 min sine signal with
amplitude greater than or equal to 1 percent should usually be recovered in
these data. Signals at 0.2--0.5 percent are not reliably recovered. Therefore,
if V0399 UMa varied at the VSX period during our observations, the modulation
was likely below about the percent level, non-sinusoidal, or hidden by the
systematic effects described above.

### 5.3 AM CVn

AM CVn gave a more promising light curve. The AAVSO sequence solution placed
the target and comparison stars directly on detected sources in the good
reference frame. The final differential light curve used 138 good target frames.
The comparison-control scatter ranged from 0.70 percent to 1.22 percent for the
AAVSO sequence stars. Most control stars were statistically consistent with
constant flux under the empirical scatter model; star 90 was the least clean
control and also saturated in nine frames.

The strongest Lomb-Scargle peak in the 5--90 min search is at 8.766 min, with
power 0.156 and false-alarm probability 0.0138. A sine fit at this period has
amplitude 0.93 percent and improves BIC by 13.47 relative to a constant model.
The wider 3--180 min search gives a consistent peak at 8.744 min with FAP
0.0140. Fixed-period fits near 17 min are weaker but still modestly favored:
17.14 min gives amplitude 0.72 percent and BIC improvement 3.86; 17.20 min
gives amplitude 0.74 percent and BIC improvement 4.36.

![Figure 2. AM CVn differential light curve. Fluxes are exposure-normalized and
divided by an AAVSO comparison-star ensemble formed after normalizing each
comparison star by its own median flux. The final two artifact frames are
excluded. The light curve shows percent-level short-timescale structure.](/Users/josephbailey/Documents/School/100/ph100-delta-scutis/analysis/project/plots/AM_CVn_sequence_X42421BZ_light_curve.png)

![Figure 3. AM CVn Lomb-Scargle periodogram from 5 to 90 min. The strongest
peak is near 8.77 min with false-alarm probability 0.014. Because this period
is close to half of the expected AM CVn timescale, it may be a harmonic or alias
of the physical orbital/superhump modulation.](/Users/josephbailey/Documents/School/100/ph100-delta-scutis/analysis/project/periodograms/AM_CVn_sequence_X42421BZ/short_5_to_90min_periodogram.png)

The interpretation is therefore stronger for AM CVn than for V0399 UMa, but it
still requires caution. The 8.76 min signal may be the first harmonic of a
roughly 17.5 min modulation rather than a distinct physical period. A longer
baseline or a second night of observations would help separate the harmonic from
the fundamental and reduce the risk of an alias.

## 6. Conclusions

Our May 20 observations produced two different outcomes.

For V0399 UMa, we do not claim a period detection. The final light curve has a
noise floor near 1.4 percent. The cataloged 15.3 min period is not recovered
significantly. Aperture-size tests, difference imaging, and two-star models all
show that the remaining V0399 residuals are strongly coupled to DN UMa or to
common-mode image systematics. The scientifically defensible conclusion is that
our data are consistent with constant flux within the systematic noise floor,
while weak sub-percent variability cannot be ruled out.

For AM CVn, we do detect a credible short-period signal. The dominant
Lomb-Scargle peak is near 8.76 min with FAP about 0.014 and amplitude about
0.93 percent. Fixed fits near 17.1 min are weaker but still moderately favored.
The most likely interpretation is that our data detect short-timescale AM CVn
variability, but the limited baseline and harmonic ambiguity prevent a precise
period determination from this dataset alone.

The main limitations of the project were calibration imperfections, tracking
drift, comparison-star saturation, and short time baseline. The most important
improvements would be to obtain a longer continuous run, use cleaner twilight or
dome flats without stars, avoid saturated comparison stars, and obtain repeated
observations on multiple nights. For V0399 UMa specifically, higher spatial
resolution or PSF photometry would be useful because DN UMa is close enough to
contaminate aperture and difference-image measurements.

## Acknowledgments

We thank the PHYSICS 100 teaching staff for telescope access, observing support,
and data-reduction guidance. We used AAVSO resources including the Variable Star
Plotter sequence for AM CVn and the VSX catalog as a reference for known
variable-star periods. TODO: add specific names and institutional acknowledgments.

## References

- Astropy Collaboration. *Lomb-Scargle Periodograms*, Astropy documentation.
  https://docs.astropy.org/en/stable/timeseries/lombscargle.html
- AAVSO. *Variable Star Plotter and comparison-star charts*.
  https://www.aavso.org/variable-star-charts
- AAVSO. *The International Variable Star Index (VSX)*.
  https://vsx.aavso.org/
- Baluev, R. V. 2008, *Assessing the statistical significance of periodogram
  peaks*, MNRAS, 385, 1279.
  https://academic.oup.com/mnras/article/385/3/1279/1010111
- Levitan, D., et al. 2013, *Orbital periods and accretion disc structure of
  four AM CVn systems*, MNRAS, 432, 2048.
  https://academic.oup.com/mnras/article/432/3/2048/1747230
- Lomb, N. R. 1976, *Least-squares frequency analysis of unequally spaced
  data*, Astrophysics and Space Science, 39, 447.
- Scargle, J. D. 1982, *Studies in astronomical time series analysis. II.
  Statistical aspects of spectral analysis of unevenly spaced data*,
  Astrophysical Journal, 263, 835.
- Solheim, J.-E. 2010, *AM CVn Stars: Status and Challenges*, PASP, 122, 1133.

## Appendix A. Individual Contributions

TODO: Fill this in before submission. The project guidelines require individual
team member contributions to be identified in an appendix.

Suggested format:

| Team member | Contributions |
|---|---|
| TODO | Observing, calibration analysis, V0399 photometry, report writing |
| TODO | AM CVn photometry, periodograms, figures, presentation preparation |
| TODO | Literature background, uncertainty analysis, editing |

## Appendix B. Reproducible Analysis Products

Key files used in this draft:

- `docs/ANALYSIS_PLAN.md`
- `docs/RESULTS.md`
- `scripts/project_pipeline.py`
- `scripts/v0399_contamination_tests.py`
- `scripts/v0399_polish_figures.py`
- `scripts/amcvn_sequence_analysis.py`
- `analysis/project/contamination_tests/V0399_UMa/plots/V0399_final_summary_panel.png`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/differential_light_curve.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/period_results.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/fixed_period_sine_summary.csv`
