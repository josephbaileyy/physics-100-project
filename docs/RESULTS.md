# V0399 UMa and AM CVn Results and Uncertainty Tests

This document summarizes the current V0399 UMa reduction results, the tests
that were run to reduce uncertainty, and the claims that are supported by the
data. It is a results companion to `docs/ANALYSIS_PLAN.md`.

## Current Best Interpretation

The May 20 V0399 data do not give a secure Lomb-Scargle period detection. The
data are useful, but the residual noise floor is still at about the `1.4%`
level for the most internally consistent combined light curve.

The most defensible combined analysis uses:

- Per-frame SEP matching for comparison stars.
- The same 8 Gaia comparison stars in both the `light_V_8s` and
  `light_V0399_V_8s` sequences.
- `flat_4s_masked_clipmean` as the lowest target-scatter flat-field candidate,
  with `no_flat` retained as an important systematic control.

With that setup:

- V0399 robust scatter: `1.385%`.
- DN/check robust scatter: `1.673%`.
- Best Lomb-Scargle peak after binning tests remains high-FAP and not secure.
- VSX period `15.3 min` is not recovered significantly in this dataset.
- Constant flux is an acceptable model when empirical check-star uncertainty is
  used with the common-comparison reduction.

## Data Used

Two V0399-related sequences were analyzed:

- `light_V_8s`: 83 frames, spanning `22.63 min`.
- `light_V0399_V_8s`: 109 frames, spanning `29.80 min`.

Together they span `54.11 min`, with a gap between the two sequences.

The `light_V_8s` sequence was confirmed to contain V0399 UMa. It required its
own source solve because the named V0399 reference coordinates were not valid
for this pointing.

## Source Solving and Photometry

Early attempts used reference-frame comparison-star coordinates shifted by the
V0399/DN pair. That worked for the first few `light_V_8s` frames but failed
after frame 5: comparison apertures landed on blank sky and the differential
flux became artificially huge.

The working solution is per-frame SEP matching:

1. Detect V0399 and DN UMa directly in every frame.
2. Use the V0399/DN pair to predict the field transform.
3. Detect sources with SEP in that same frame.
4. Match each comparison star to the nearest SEP source.
5. Perform aperture photometry at the per-frame matched positions.

This fixed the main solve failure:

- `light_V_8s`: `83/83` good frames with `7-8` comparison stars.
- `light_V0399_V_8s`: `107/109` good frames; 2 rejected because V0399 saturated.

Useful per-sequence products:

- `analysis/project/photometry/V0399_UMa_V8s_per_frame_sep/V8s_per_frame_sep_light_curve.csv`
- `analysis/project/photometry/V0399_UMa_named_per_frame_sep/V0399_named_per_frame_sep_light_curve.csv`

The old direct-tracking plot and the combined per-frame SEP light curve were
removed from the shared products because they are easy to misuse: direct
tracking failed on later `V_8` frames, and the combined per-frame SEP product
used different comparison ensembles for the two sequences.

## Comparison-Star Selection

Several comparison ensembles were tested. Optimizing each sequence separately
reduced scatter, but it introduced a suspicious bridge/trend between the two
segments because the two sequences used different comparison ensembles.

The most consistent combined setup uses all 8 Gaia comparison stars common to
both sequences:

- `gaia_774271076054448768_src2`
- `gaia_774272862760870784_src1`
- `gaia_774287774887292672_src5`
- `gaia_774288805679445888_src7`
- `gaia_774289041901756032_src4`
- `gaia_774289488578357504_src6`
- `gaia_774290626745586944_src8`
- `gaia_774384016513181696_src3`

Shared comparison-set tests:

| Comparison set | V0399 named scatter | V_8 scatter | named/V_8 median ratio | trend p-value |
|---|---:|---:|---:|---:|
| all common 8 | `1.37%` | `1.58%` | `1.003` | `0.51` |
| `src2,src1,src5,src4` | `1.09%` | `1.37%` | `1.001` | `0.020` |
| `src2,src7,src8` | `2.16%` | `1.76%` | `0.998` | `0.53` |
| bright common 5 | `1.69%` | `2.39%` | `1.008` | `0.00025` |

Conclusion: the all-common-8 set is preferred for combined claims because it
removes the artificial inter-segment trend. Smaller optimized sets can reduce
scatter, but they are more likely to induce segment-dependent behavior.

Relevant products:

- `analysis/project/photometry/V0399_UMa_shared_comparison_sets/shared_comparison_set_summary.csv`
- `analysis/project/plots/V0399_UMa_shared_comparison_set_tests.png`

## Calibration Choices Tested

The following flat-field candidates were tested with the same common-8
comparison ensemble:

- `flat_0p3_late_bias`: dome-light `0.3s` flats minus late clean bias.
- `flat_0p3_early_bias`: dome-light `0.3s` flats minus early dome-light bias.
- `flat_4s_masked`: star-masked `4s` flats, median combine.
- `flat_6s_masked`: star-masked `6s` flats, median combine.
- `no_flat`: diagnostic reduction with no flat correction.
- `flat_4s_masked_clipmean`: star-masked `4s` flats, sigma-clipped mean.
- `flat_6s_masked_clipmean`: star-masked `6s` flats, sigma-clipped mean.

The masking used for star-contaminated flats:

1. Bias/dark subtract each flat.
2. Normalize each flat by its median.
3. Median-filter to estimate the broad illumination pattern.
4. Subtract the smooth model to isolate compact residuals.
5. Mask bright residuals above a sigma threshold.
6. Dilate the mask to include star wings.
7. Combine the masked flats.

The original `flat_4s_masked` and `flat_6s_masked` used a masked median. The
new `*_clipmean` tests used a masked, sigma-clipped mean.

Common-8 photometry by calibration:

| Calibration | Combined V0399 scatter | Combined DN/check scatter |
|---|---:|---:|
| `flat_4s_masked_clipmean` | `1.385%` | `1.673%` |
| `flat_4s_masked` | `1.392%` | `1.732%` |
| `no_flat` | `1.394%` | `1.574%` |
| `flat_0p3_late_bias` | `1.424%` | `1.757%` |
| `flat_6s_masked_clipmean` | `1.425%` | `1.599%` |
| `flat_0p3_early_bias` | `1.429%` | `1.744%` |
| `flat_6s_masked` | `1.444%` | `1.654%` |

Conclusions:

- The `4s` masked sigma-clipped mean flat is currently the best target-scatter
  calibration.
- The improvement over the median `4s` flat is small: `1.392%` to `1.385%`.
- `no_flat` is surprisingly competitive and should remain a systematic check.
- The `0.3s` early-bias and late-bias reductions are essentially identical.
- Calibration choice changes the answer, but not enough to reveal the VSX
  period.

Relevant products:

- `analysis/project/calibration/flat_candidate_photometry_qc/flat_candidate_common8_photometry_summary_with_clipmean.csv`
- `analysis/project/plots/V0399_flat_candidate_clipmean_comparison.png`

## Calibration Uncertainty Quantification

The frame-by-frame spread in V0399 normalized differential flux across the real
flat candidates is:

| Sequence | Median calibration sigma | 90th percentile sigma | Median calibration range |
|---|---:|---:|---:|
| V_8 generic | `0.134%` | `0.403%` | `0.304%` |
| named V0399 | `0.117%` | `0.311%` | `0.252%` |
| combined | `0.119%` | `0.353%` | `0.257%` |

Including `no_flat` increases the combined median calibration sigma to
`0.170%`.

The dome-light-bias test is negligible:

| Metric | Combined value |
|---|---:|
| median `0.3s early-bias - late-bias` effect | `0.001%` |
| robust scatter of early/late difference | `0.003%` |
| 90th percentile absolute difference | `0.0065%` |

Flat pixel-level uncertainty across real flat candidates:

- median per-pixel flat sigma: `0.259%`
- 90th percentile: `0.477%`
- 99th percentile: `0.687%`

Conclusion: calibration uncertainty is real but secondary. It contributes
roughly `0.1-0.4%` pointwise, compared with a final light-curve scatter near
`1.4-1.7%`.

Relevant products:

- `analysis/project/calibration/flat_candidate_photometry_qc/calibration_pointwise_uncertainty_summary.csv`
- `analysis/project/calibration/flat_candidate_photometry_qc/calibration_bias_choice_uncertainty_summary.csv`
- `analysis/project/calibration/flat_candidate_photometry_qc/flat_pixel_uncertainty_summary.csv`

## Tracking and Position-Dependent Systematics

Tracking drift was tested by correlating residuals with measured target x/y
position, time, ensemble flux, and solve-quality metrics.

Key correlations:

- V_8 residual vs time: Pearson `r = 0.276`, p `0.012`.
- V_8 residual vs x drift: `r = 0.258`, p `0.019`.
- V_8 residual vs y drift: `r = -0.217`, p `0.048`.
- Named V0399 residual vs x/y/time: not significant, p about `0.18-0.23`.
- Combined residual vs x drift: weak, p `0.010`.

Linear regression results:

- Position-only explains about `5%` in-sample variance.
- Leave-one-out predictive power for position-only is approximately zero or
  negative.
- Time-only is similarly weak.

Conclusion: tracking/calibration interaction may contribute to the noise floor,
especially in V_8, but measured x/y drift does not explain most of the residual
scatter.

Relevant products:

- `analysis/project/systematics/V0399_tracking_calibration/residual_feature_correlations.csv`
- `analysis/project/systematics/V0399_tracking_calibration/residual_regression_model_comparison.csv`
- `analysis/project/systematics/V0399_tracking_calibration/residuals_vs_tracking_calibration.png`

## Constant-Flux Tests

Using the separately optimized balanced subset, V0399 was formally
nonconstant, but DN/check also showed nonconstant behavior, which weakens the
claim:

| Series | N | empirical sigma | reduced chi2 | p-value constant |
|---|---:|---:|---:|---:|
| V0399 combined | 190 | `1.19%` | `1.41` | `1.7e-4` |
| DN/check combined | 142 | `1.19%` | `1.35` | `3.9e-3` |

With the all-common-8 comparison set, the constant model is acceptable under
the empirical DN/check noise model. This is another reason to prefer the
common-comparison result for claims.

HD 103187 was also tested as a target-like control. With a balanced comparison
ensemble excluding HD 103187:

- HD 103187 scatter: `2.07%`
- reduced chi2: `1.16`
- p-value for constant flux: `0.125`

Conclusion: HD 103187 is consistent with constant flux in the safer test, while
V0399 has at most marginal extra structure.

Relevant products:

- `analysis/project/photometry/V0399_UMa_balanced_comparison_subset_chi2.csv`
- `analysis/project/photometry/HD_103187_constant_flux_chi2.csv`

## Period and Variable-Star Model Tests

### Lomb-Scargle

The optimized/balanced comparison light curve gave:

- Best short-period peak: `4.68 min`
- amplitude: `0.40%`
- RMS: `1.39%`
- FAP: `1.0`

A Delta-Scuti-ish range search gave:

- Best period: `29.52 min`
- amplitude: `0.32%`
- RMS: `1.40%`
- FAP: `1.0`

Conclusion: no secure Lomb-Scargle detection.

### VSX Period

VSX gives V0399 UMa period `15.3 min`. A fixed-period fit to the balanced
subset gave:

| Model | Period | amplitude | RMS | reduced chi2 | improvement vs constant |
|---|---:|---:|---:|---:|---:|
| constant | | `0%` | `1.409%` | `1.41` | |
| VSX fixed period | `15.3 min` | `0.284%` | `1.395%` | `1.40` | p `0.069` |
| best LS peak | `4.68 min` | `0.401%` | `1.381%` | `1.37` | p `0.0048` |

The VSX period is only weakly hinted, not recovered.

### Other Variable-Star Models

On the all-common-8 light curve, these models were compared:

- constant per sequence
- sine waves
- two-harmonic sine waves
- broad eclipse/box-like model
- linear trend

The constant model is already acceptable with empirical uncertainties:

- constant per sequence BIC: `145.5`
- broad box model BIC: `145.0`
- best sine at `5.41 min` BIC: `149.2`
- sine at `15.3 min` BIC: `155.6`

The broad box model is not physically compelling: it has period `3.105 min`
and duration `2.40 min`, so it is better viewed as a flexible residual-shape
fit than as a real eclipsing model.

Conclusion: no tested variable-star model is convincing. The best scientific
statement is that V0399 is consistent with constant flux within the current
systematic noise floor, with possible weak sub-percent structure that this
dataset cannot confirm.

Relevant products:

- `analysis/project/periodograms/V0399_UMa_balanced_comparison_subset/vsx_15p3min_fixed_period_fit_summary.csv`
- `analysis/project/model_fits/V0399_UMa_shared_common8_variable_models/variable_model_comparison_common8.csv`

## Binning / Combining Frames

Adjacent-frame binning was tested on the current best calibration/common-8
light curve.

| Binning | V0399 scatter | DN/check scatter | Best LS period | LS FAP |
|---|---:|---:|---:|---:|
| unbinned | `1.38%` | `1.67%` | `5.41 min` | `0.99995` |
| 2-frame | `1.09%` | `1.26%` | `5.43 min` | `0.9566` |
| 3-frame | `0.87%` | `1.25%` | `5.43 min` | `0.9999` |
| 5-frame | `0.57%` | `0.94%` | `5.29 min` | `0.882` |
| 8-frame | `0.66%` | `0.87%` | `3.01 min` | `0.279` |

Conclusion: binning makes the curve cleaner visually, but it does not produce a
credible Lomb-Scargle detection. It is useful for plots, not for a stronger
period claim.

Relevant products:

- `analysis/project/model_fits/V0399_UMa_binning_tests/binning_model_summary.csv`
- `analysis/project/model_fits/V0399_UMa_binning_tests/binning_light_curves.png`

## Blending, Difference Imaging, and Injection-Recovery

The final robustness suite uses the current preferred reduction:
`flat_4s_masked_clipmean`, the all-common-8 comparison ensemble, and both
`light_V_8s` plus `light_V0399_V_8s`.

### Aperture / Annulus Blending Test

Photometry was rerun over aperture radii `3, 4, 5, 6, 7, 8, 10 px` and sky
annuli `8-14`, `10-18`, `14-24`, and `18-30 px`.

Median behavior by aperture radius:

| Radius | V0399 scatter | DN scatter | HD 103187 scatter | V0399-DN corr. | 15.3 min amp. | LS FAP |
|---:|---:|---:|---:|---:|---:|---:|
| 3 px | `3.65%` | `3.95%` | `5.38%` | `0.17` | `0.326%` | `0.99998` |
| 5 px | `1.39%` | `1.59%` | `2.06%` | `0.26` | `0.137%` | `1.00000` |
| 8 px | `0.67%` | `0.80%` | `0.90%` | `0.66` | `0.121%` | `1.00000` |
| 10 px | `0.68%` | `0.63%` | `0.71%` | `0.78` | `0.048%` | `0.99999` |

The lowest scatter occurs for larger apertures, but those same apertures
greatly increase the V0399-DN correlation. That is a warning sign for
common-mode residuals, blended wings, imperfect background treatment, or
registration/flat-field systematics. It does not support a clean independent
V0399 signal.

### Difference-Image Photometry

A separate reference image was built for each pointing, frames were shifted
onto the reference using the per-frame V0399/DN solve, median background/flux
scaling was applied, and residual flux was measured at V0399 and DN with
`r = 4, 5, 6 px`.

The V0399 and DN difference-image residuals are extremely correlated:

| Residual aperture | N | V0399-DN residual corr. |
|---:|---:|---:|
| 4 px | 192 | `0.917` |
| 5 px | 192 | `0.916` |
| 6 px | 192 | `0.916` |

This is not the expected signature of a clean residual at only V0399. It is
more consistent with shared subtraction residuals or nearby-pair contamination.

### Two-Star Contamination Model

The two-star mixture model compared constant V0399, V0399 mixed with DN, fixed
period components, and fixed-period components plus DN. The best BIC model is
`5.41 min sine + DN`, but it improves only slightly over `V0399 mixed with DN`:

| Model | RMS | BIC | DN leakage beta |
|---|---:|---:|---:|
| `5.41 min sine + DN` | `1.372%` | `-1193.2` | `0.229` |
| `V0399 mixed with DN` | `1.425%` | `-1192.3` | `0.260` |
| constant V0399 | `1.481%` | `-1186.4` | |
| `15.3 min sine + DN` | `1.424%` | `-1182.7` | `0.257` |
| `15.3 min sine` | `1.478%` | `-1177.1` | |

Conclusion: DN/common-mode terms explain as much or more than a physically
interesting fixed `15.3 min` component.

### Fixed-Period and Non-Sinusoidal Tests

At the VSX period, V0399 is not preferred over constant:

- Best V0399 fixed-period model overall: `5.41 min` sawtooth, BIC improvement
  `-5.82`, amplitude proxy `1.25%`.
- V0399 `15.3 min` sawtooth: BIC penalty `+4.29`, amplitude proxy `0.38%`.
- V0399 `15.3 min` Gaussian dip: BIC penalty `+13.09`, amplitude proxy
  `0.61%`.
- HD 103187 and DN do not show a credible `15.3 min` preference either.

The flexible/non-sinusoidal fits can describe small residual structure, but the
only V0399 BIC-favored shape is a `5.41 min` sawtooth. Because the broad
Lomb-Scargle FAP remains near 1 and the aperture/difference-image tests show
strong common-mode behavior, this is best treated as descriptive residual
fitting, not a variable-star classification.

### Injection-Recovery

Synthetic `15.3 min` signals were injected into V0399-like residuals, HD
103187, and five artificial/control light curves using the real cadence.

Recovery summary:

| Waveform | Amplitude/depth | Fixed-period BIC recoveries | LS-near-15.3 recoveries |
|---|---:|---:|---:|
| sine | `0.2%` | `0/7` | `0/7` |
| sine | `0.5%` | `3/7` | `1/7` |
| sine | `1.0%` | `7/7` | `5/7` |
| sine | `2.0%` | `7/7` | `7/7` |
| Gaussian dip | `0.5%` | `1/7` | `0/7` |
| Gaussian dip | `1.0%` | `1/7` | `1/7` |
| Gaussian dip | `2.0%` | `7/7` | `2/7` |
| sawtooth | `0.5%` | `3/7` | `0/7` |
| sawtooth | `1.0%` | `6/7` | `3/7` |
| sawtooth | `2.0%` | `7/7` | `6/7` |

This means the dataset should recover a coherent `15.3 min` sine at about
`1%` amplitude and very strongly at `2%`, but it cannot reliably rule out
`0.2-0.5%` sinusoidal variability or shallow non-sinusoidal/eclipselike
signals.

Relevant products:

- `analysis/project/contamination_tests/V0399_UMa/aperture_annulus_grid_summary.csv`
- `analysis/project/contamination_tests/V0399_UMa/difference_image_residual_light_curves.csv`
- `analysis/project/contamination_tests/V0399_UMa/difference_image_residual_correlation.csv`
- `analysis/project/contamination_tests/V0399_UMa/two_star_contamination_model_summary.csv`
- `analysis/project/contamination_tests/V0399_UMa/injection_recovery_summary.csv`
- `analysis/project/contamination_tests/V0399_UMa/fixed_period_model_summary.csv`
- `analysis/project/contamination_tests/V0399_UMa/non_sinusoidal_model_summary.csv`
- `analysis/project/contamination_tests/V0399_UMa/plots/aperture_radius_dependence.png`
- `analysis/project/contamination_tests/V0399_UMa/plots/v0399_dn_correlation_vs_aperture.png`
- `analysis/project/contamination_tests/V0399_UMa/plots/difference_image_contact_sheet.png`
- `analysis/project/contamination_tests/V0399_UMa/plots/fixed_15p3min_folded.png`
- `analysis/project/contamination_tests/V0399_UMa/plots/injection_recovery_heatmap.png`
- `analysis/project/contamination_tests/V0399_UMa/plots/model_comparison_bic.png`

Final contamination/detectability conclusion: the May 20 data do not support a
robust intrinsic V0399 period detection. The signal is aperture-sensitive,
strongly correlated with DN/common residuals in difference imaging, and not
improved by the fixed `15.3 min` period. A `>=1%` coherent sinusoid should often
be recoverable, so if the VSX-period modulation was present on May 20 it was
likely below about the percent level, non-sinusoidal/shallow, or hidden by
systematics from blending, subtraction, calibration, and tracking.

## What Reduced the Uncertainty

Most helpful:

- Per-frame SEP source matching instead of shifted reference apertures.
- Using the same common comparison stars across the two pointings for combined
  claims.
- Star-masked `4s` flats, especially the sigma-clipped mean version.
- Keeping DN UMa as a check star and geometry anchor.

Somewhat helpful:

- Curated comparison subsets for per-sequence diagnostics.
- Binning adjacent frames for visualization.

Not helpful or not worth using as the main conclusion:

- Using early dome-light bias with the `0.3s` dome-light flats.
- Letting each sequence use its independently optimized comparison ensemble for
  combined physical claims.
- Treating the linear bridge between the two segments as real stellar behavior.
- Claiming the `15.3 min` VSX period from this dataset alone.
- Treating aperture-dependent or DN-correlated residual structure as intrinsic
  V0399 variability.

## Current Recommended V0399 Products

Use these for the final V0399 writeup:

- Calibration comparison:
  `analysis/project/calibration/flat_candidate_photometry_qc/flat_candidate_common8_photometry_summary_with_clipmean.csv`
- Best common-comparison calibration product:
  `analysis/project/calibration/flat_candidate_photometry_qc/flat_4s_masked_clipmean_common8_light_curve.csv`
- Shared comparison set diagnostic:
  `analysis/project/photometry/V0399_UMa_shared_comparison_sets/shared_comparison_set_summary.csv`
- Binning diagnostic:
  `analysis/project/model_fits/V0399_UMa_binning_tests/binning_model_summary.csv`
- Variable model comparison:
  `analysis/project/model_fits/V0399_UMa_shared_common8_variable_models/variable_model_comparison_common8.csv`
- Contamination/detectability suite:
  `analysis/project/contamination_tests/V0399_UMa/`
- Final summary panel:
  `analysis/project/contamination_tests/V0399_UMa/plots/V0399_final_summary_panel.png`

## AM CVn Sequence Analysis

AM CVn has now been reduced with the same preferred May 20 calibration,
`flat_4s_masked_clipmean`, and the AAVSO sequence chart `X42421BZ`.

Important implementation notes:

- The AAVSO VSP JSON was saved locally as
  `analysis/project/aavso/X42421BZ.json`.
- The chart contains sequence labels `90, 96, 112, 125, 143, 144`; it does not
  include a `99` label in the downloaded chart response.
- The generic AM CVn reference-frame selection initially picked one of the bad
  final low-count artifact frames. That was replaced with a manually selected
  good mid-run reference frame near `59.6 min`.
- The final two AM CVn frames are flagged as `frame_artifact` and excluded.
- Star `90` saturated in `9` frames and is excluded only for those frames.
- Differential photometry uses exposure-normalized flux rates, then normalizes
  each comparison star by its own median before forming the ensemble. This is
  necessary because the AAVSO sequence spans `V=8.97` to `V=14.36`.

Sequence pixel solution:

| Label | Role | V mag | x | y |
|---|---|---:|---:|---:|
| target | target | | `502.9` | `1283.9` |
| 90 | comparison | `8.97` | `732.3` | `201.4` |
| 96 | comparison | `9.551` | `1808.5` | `1067.7` |
| 112 | comparison | `11.235` | `1273.2` | `1329.0` |
| 125 | comparison | `12.454` | `389.7` | `1174.7` |
| 143 | comparison | `14.29` | `826.9` | `1385.3` |
| 144 | comparison | `14.357` | `285.1` | `1362.1` |

The best tested circular aperture for AM CVn was `r=6 px` with sky annulus
`12-20 px`.

AM CVn comparison-control scatter:

| Control | N | robust scatter | constant p-value |
|---|---:|---:|---:|
| 90 | 129 | `0.70%` | `0.044` |
| 96 | 138 | `0.81%` | `0.46` |
| 125 | 138 | `0.84%` | `0.53` |
| 112 | 138 | `0.96%` | `0.99` |
| 143 | 138 | `1.08%` | `0.23` |
| 144 | 138 | `1.22%` | `0.99` |

AM CVn period search:

| Run | N | best period | LS FAP | sine amp | sine BIC improvement |
|---|---:|---:|---:|---:|---:|
| `5-90 min` | 138 | `8.766 min` | `0.0138` | `0.93%` | `-13.47` |
| `3-180 min` | 138 | `8.744 min` | `0.0140` | `0.93%` | `-13.43` |
| fixed `17.14 min` | 138 | `17.14 min` | | `0.72%` | `-3.86` |
| fixed `17.20 min` | 138 | `17.20 min` | | `0.74%` | `-4.36` |

Interpretation: AM CVn shows a much more credible short-period signal than
V0399 in these data. The strongest Lomb-Scargle peak is near `8.75-8.77 min`,
while a fixed `17.1-17.2 min` sine is weaker but still modestly favored over a
constant model. Because `8.76 min` is close to half of `17.5 min`, this should
be treated as a possible harmonic/alias issue rather than immediately as a
distinct physical period.

### AM CVn period uncertainty and literature comparison

The formal target aperture errors underestimate the full point-to-point scatter.
The median formal normalized error is `0.90%`, while the robust target scatter is
`1.51%`. Propagating an estimated comparison-ensemble term raises the median
error only to about `0.98-1.00%`, so the remaining high reduced chi2 is best
explained by a combination of intrinsic AM CVn flickering, changing double-wave
shape, and unmodeled systematic scatter rather than a clearly incorrect
differential-photometry reduction.

Model tests on the 138 retained AM CVn points:

| Model | chi2 | reduced chi2 | BIC improvement vs constant | RMS |
|---|---:|---:|---:|---:|
| constant flux | `457.85` | `3.34` | `0.00` | `1.66%` |
| sine at LS best, `8.766 min` | `383.90` | `2.84` | `-64.10` | `1.52%` |
| sine at `17.14 min` | `412.62` | `3.06` | `-35.37` | `1.58%` |
| double wave at `17.14 min` | `359.08` | `2.70` | `-79.06` | `1.48%` |

The F-tests and BIC favor a periodic model over constant flux, and favor a
double-wave harmonic model over a single `17.14 min` sine. However, all reduced
chi2 values remain high, so the simple models do not explain all of the observed
scatter.

Allowing the period to vary gives:

| Fit | Best period | formal delta-chi2 interval | scaled-error interval |
|---|---:|---:|---:|
| sine near LS peak | `8.750 min` | `-0.042/+0.043 min` | `+/-0.072 min` |
| double-wave fundamental | `17.586 min` | `-0.089/+0.090 min` | `-0.146/+0.150 min` |

Bootstrap Lomb-Scargle tests are bimodal between the `~8.8 min` harmonic and the
`~17-19 min` family. Row-resampling gives a median FAP `4.26e-4`, with `94.7%`
of trials below `0.05`; perturbing the data by the empirically scaled errors
gives a median FAP `0.266`, with only `21.4%` below `0.05`. Thus the data
support AM CVn variability, but the exact period assignment is not robust at the
precision needed for a period-derivative claim.

Conditioning the bootstrap samples on which period family is selected gives
tighter within-family spreads, but these should not be confused with the global
period uncertainty because they condition away the alias ambiguity:

| Bootstrap | Family | Weight | conditional period | conditional FAP |
|---|---:|---:|---:|---:|
| row resample | half-period | `0.608` | `8.766 -0.085/+0.070 min` | median `3.35e-4` |
| row resample | doubled half-period | `0.608` | `17.532 -0.171/+0.139 min` | median `3.35e-4` |
| row resample | long-period family | `0.373` | `18.289 -0.477/+0.678 min` | median `6.38e-4` |
| scaled perturbation | half-period | `0.556` | `8.766 -0.069/+0.070 min` | median `0.228` |
| scaled perturbation | doubled half-period | `0.556` | `17.532 -0.137/+0.139 min` | median `0.228` |
| scaled perturbation | long-period family | `0.433` | `18.424 -0.591/+0.631 min` | median `0.317` |

A two-Gaussian mixture fit to the bootstrap period distribution was also tested.
It gives a compact half-period component near `8.75 min`, but the long-period
component is broad and non-Gaussian because it contains aliases/tails. Therefore
the family-split percentile intervals above are preferred for reporting
conditional bootstrap widths.

Patterson et al. (1979) measured the 1978 photometric/superhump-like period as
`1051.212 +/- 0.015 s` and reported `dP/dt = +0.012 +/- 0.003 s/yr`. Projected
from their 1978 epoch to the midpoint of these observations
(`JD 2461181.77`, 2026-05-21 UTC), that predicts only
`1051.790 +/- 0.145 s`, a change of `0.578 +/- 0.145 s`. Our best double-wave
fundamental is `17.586 min = 1055.18 s`, with a scaled uncertainty of about
`8.9 s`. This is consistent with both Patterson's 1978 photometric period and
the period predicted by their quoted derivative, but the uncertainty is roughly
`15x` larger than the expected secular change. Therefore these data cannot
confirm or refute the Patterson period increase.

Roelofs et al. (2006) distinguish the true spectroscopic orbital period
`1028.7322 +/- 0.0003 s` from the main `~1051 s` photometric superhump. Our
photometric double-wave period is consistent with the superhump family, but it
is not a precision orbital-period measurement and should not be used as a test
of the spectroscopic orbital period.

Kruszewski & Semeniuk (1992) reanalyzed Smak's 1962 observations and likewise
separated the dominant `1051 s` family, mostly seen through the `525 s` harmonic,
from another weaker `1023` or `1011 s` mode. They concluded that the period of
the dominant `525 s` harmonic could not be uniquely determined from the 1962
data, probably because of slow changes over the four-month observing interval.
They preferred `1023.02 s` for the other mode, while noting that Patterson et al.
(1992) obtained `1023.44 s` and that the difference could be a cycle-count alias.

Our LS half-period is `524.985 s`; doubled, it is `1049.970 s` with a scaled
uncertainty of `8.64 s`. The direct double-wave fundamental fit gives
`1055.179 s` with scaled uncertainty `8.89 s`. Both are consistent with the
`~1051 s` superhump family discussed by Kruszewski & Semeniuk. They are
inconsistent with `1023.02 s` at about `3.1-3.6 sigma` using the scaled
uncertainties, but that is not evidence for a secular period change because
`1023 s` is a different mode from the dominant `1051 s`/`525 s` superhump
component measured here.

### Exploratory AM CVn period-modeling pass

A more complete exploratory pass was run with data-driven fits only, followed by
literature comparison. The baseline is `111.20 min` with 138 retained points and
median cadence `48.55 s`. This covers about `12.7` cycles of the `525 s`
harmonic, `6.3` cycles of the `1051 s` superhump, and only `0.139` cycles of the
`13.38 h` disk precession.

The best BIC model is a free-period double-wave model:

| Model | Period | chi2 | reduced chi2 | BIC improvement vs constant | RMS |
|---|---:|---:|---:|---:|---:|
| double wave, free period | `17.5866 min` | `336.46` | `2.53` | `-101.68` | `1.43%` |
| sine, LS half-period | `8.7422 min` | `383.79` | `2.84` | `-64.20` | `1.53%` |
| sine, `1051.2 s` superhump | `17.5200 min` | `403.99` | `2.99` | `-44.00` | `1.56%` |
| sine, `1028.7322 s` orbital | `17.1455 min` | `412.49` | `3.06` | `-35.50` | `1.58%` |

The double-wave model strongly improves the fit, but reduced chi2 remains above
1, so extra flickering/waveform variation is still present.

Windowed Lomb-Scargle scans find:

| Window | Best period | FAP | Interpretation |
|---|---:|---:|---|
| `5-12 min` | `8.744 min` | `0.0085` | robust half-period harmonic |
| `7-10 min` | `8.742 min` | `0.0086` | robust half-period harmonic |
| `14-22 min` | `18.366 min` | `0.043` | weaker long-period family |
| `3-180 min` | `8.744 min` | `0.0084` | global best is harmonic |
| `115-120 s` | `119.6 s` | `1.0` | not detected |
| `26.3 s` window | `23.8 s` | `0.296` | not detected; below cadence/exposure limit |

Residual periodograms after subtracting the best double-wave model do not show a
credible remaining signal: the `115-120 s` residual search has FAP `1.0`, and
the `26.3 s` residual search has FAP `0.997`.

Ten-thousand-trial bootstrap tests show that the period estimate remains
bimodal:

| Bootstrap | FAP median | FAP < 0.05 | half-family fraction | long-family fraction |
|---|---:|---:|---:|---:|
| row resample | `4.64e-4` | `0.950` | `0.596` | `0.387` |
| perturb formal errors | `0.054` | `0.482` | `0.646` | `0.353` |
| perturb scaled errors | `0.253` | `0.208` | `0.573` | `0.411` |
| 10 min block resample | `1.02e-4` | `0.973` | `0.613` | `0.386` |

Conditional family intervals are tighter but do not remove the global alias
ambiguity. For the scaled-error perturbation bootstrap, the half-period family
is `8.766 -0.102/+0.070 min`; doubled, this is
`17.532 -0.205/+0.139 min`. The long-period family is
`18.424 -0.591/+0.631 min`.

The detectability limits are important:

- `26.3 s`: not reliable; shorter than the exposures and below the cadence
  Nyquist, with median exposure-smearing attenuation only `0.21`.
- `115-120 s`: marginally sampled at only `2.4` samples per cycle and not
  detected in the residual periodogram.
- `13.38 h` precession: not measurable from this run because less than one
  fifth of a cycle is covered.

Final AM CVn interpretation from the exploratory pass: the data favor a
photometric superhump-family signal dominated by the `~525 s` harmonic /
`~1051 s` double-wave family. The spectroscopic `1028.7 s` orbital period is not
the dominant photometric signal in these data, but this does not contradict the
literature because the orbital period is spectroscopic and the photometry is
expected to be superhump-dominated. The present data are not precise enough to
test the Patterson period derivative or to measure disk precession.

Relevant AM CVn products:

- `analysis/project/aavso/X42421BZ.json`
- `analysis/project/aavso/X42421BZ_sequence_pixel_table.csv`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/differential_light_curve.csv`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/aperture_photometry.csv`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/comparison_control_summary.csv`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/photometry_flag_summary.csv`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/AM_CVn_sequence_reference_overlay.png`
- `analysis/project/photometry/AM_CVn_sequence_X42421BZ/AM_CVn_early_late_frame_contact_sheet.png`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/period_results.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/fixed_period_sine_summary.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/model_chi2_f_bic_summary_refined.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/period_delta_chi2_intervals.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/bootstrap_period_summary.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/bootstrap_period_family_conditional_intervals.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/bootstrap_period_family_sklearn_gmm_summary.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/literature_period_change_test_summary.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/literature_period_comparisons.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/literature_1992_kruszewski_comparison.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/exploratory_model_comparison.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/exploratory_lomb_scargle_window_summary.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/exploratory_bootstrap_period_summary.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/exploratory_bootstrap_family_intervals.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/exploratory_residual_periodogram_summary.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/exploratory_signal_detectability_summary.csv`
- `analysis/project/periodograms/AM_CVn_sequence_X42421BZ/exploratory_literature_comparison.csv`
- `analysis/project/plots/AM_CVn_sequence_X42421BZ_light_curve.png`
- `analysis/project/plots/AM_CVn_period_uncertainty_diagnostics.png`
- `analysis/project/plots/AM_CVn_exploratory_periodogram_windows.png`
- `analysis/project/plots/AM_CVn_exploratory_folded_models.png`
- `analysis/project/plots/AM_CVn_exploratory_bootstrap_family_histograms.png`
- `analysis/project/plots/AM_CVn_exploratory_residual_periodograms.png`

## Final Suggested Claim

The May 20 V0399 UMa observations support a careful null or marginal-variability
claim, not a period detection:

> We reduced the May 20 V0399 UMa frames with per-frame SEP source matching and
> tested multiple calibration and comparison-star configurations. The lowest
> target-scatter reduction used a star-masked sigma-clipped mean `4s` flat, but
> calibration choice contributed only about `0.1-0.4%` pointwise uncertainty.
> Using the same 8 Gaia comparison stars across both pointings gave the most
> internally consistent combined light curve, with V0399 scatter near `1.4%`.
> Lomb-Scargle, fixed-period VSX fitting at `15.3 min`, binned light curves, and
> alternate variable-star models did not produce a secure period detection.
> Aperture tests and difference-image residuals show strong coupling between
> V0399 and nearby DN/common-mode structure, so the remaining weak signals are
> not robustly intrinsic. Injection-recovery shows that coherent `>=1%`
> sinusoidal modulation at `15.3 min` should usually be detectable, while
> `0.2-0.5%` or shallow non-sinusoidal modulation could be missed. The data are
> consistent with constant flux within the empirical systematic noise floor,
> though weak low-amplitude variability cannot be ruled out.
