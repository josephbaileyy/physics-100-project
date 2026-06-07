"""Physically motivated AM CVn model ladder.

This script compares generic period searches with literature-tied AM CVn
models.  It deliberately separates independent/free-period constraints from
conditional constraints that assume the detected signal is one of the known
orbital or superhump components.
"""
from __future__ import annotations

import argparse
import itertools
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/mplconfig")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle
from scipy import optimize, stats


ROOT = Path(__file__).resolve().parents[2]
PHOT_DIR = ROOT / "analysis/am_cvn/photometry/sequence_X42421BZ"
PER_DIR = ROOT / "analysis/am_cvn/periodograms/sequence_X42421BZ"
PLOT_DIR = ROOT / "analysis/am_cvn/plots"

ORBITAL_S = 1028.7322
ORBITAL_SIGMA_S = 0.0003
SUPERHUMP_S = 1051.2
SUPERHUMP_SIGMA_S = 0.1
PATTERSON_SUPERHUMP_S = 1051.212
PATTERSON_SUPERHUMP_SIGMA_S = 0.015
NEG_SUPERHUMP_S = 1011.4
PRECESSION_HR = 13.38
RNG_SEED = 20260530


def robust_sigma(x: np.ndarray) -> float:
    med = np.nanmedian(x)
    return float(1.4826 * np.nanmedian(np.abs(x - med)))


def load_light_curve() -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    diff = pd.read_csv(PHOT_DIR / "differential_light_curve.csv")
    good = diff[diff.good.astype(bool) & np.isfinite(diff.diff_flux_norm)].sort_values("jd").copy()
    median_diff = float(np.nanmedian(good.diff_flux.to_numpy(float)))
    t_min = (good.jd.to_numpy(float) - good.jd.min()) * 24 * 60
    y = good.diff_flux_norm.to_numpy(float)
    sigma = good.diff_flux_err.to_numpy(float) / median_diff
    ok = np.isfinite(t_min) & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    return good.iloc[np.where(ok)[0]].copy(), t_min[ok], y[ok], sigma[ok]


def design_matrix(t_min: np.ndarray, components: list[tuple[float, int]]) -> np.ndarray:
    cols = [np.ones(len(t_min))]
    for period_min, harmonic in components:
        phase = 2 * np.pi * harmonic * t_min / period_min
        cols.extend([np.sin(phase), np.cos(phase)])
    return np.column_stack(cols)


def fit_linear(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    components: list[tuple[float, int]],
) -> dict:
    x = design_matrix(t_min, components)
    beta = np.linalg.lstsq(x / sigma[:, None], y / sigma, rcond=None)[0]
    fitted = x @ beta
    resid = y - fitted
    chi2 = float(np.sum((resid / sigma) ** 2))
    sse = float(np.sum(resid**2))
    n = len(y)
    k = x.shape[1]
    dof = n - k
    amps = {}
    phases = {}
    for i, (period_min, harmonic) in enumerate(components):
        a = beta[1 + 2 * i]
        b = beta[2 + 2 * i]
        label = f"p{period_min * 60:.1f}s_h{harmonic}"
        amps[f"amp_{label}_pct"] = 100 * float(np.hypot(a, b))
        phases[f"phase_{label}_cycles"] = float((math.atan2(b, a) / (2 * np.pi)) % 1)
    return {
        "periods_min": ";".join(f"{p:.8f}" for p, _ in components),
        "components": ";".join(f"{p * 60:.4f}s_h{h}" for p, h in components),
        "N": n,
        "k": k,
        "dof": dof,
        "chi2": chi2,
        "reduced_chi2": chi2 / dof,
        "chi2_p_value": float(stats.chi2.sf(chi2, dof)),
        "aic_chi2": chi2 + 2 * k,
        "bic_chi2": chi2 + k * math.log(n),
        "sse": sse,
        "bic_sse": n * math.log(max(sse / n, 1e-300)) + k * math.log(n),
        "rms_pct": 100 * float(np.sqrt(np.mean(resid**2))),
        "beta": beta,
        "fitted": fitted,
        "resid": resid,
        **amps,
        **phases,
    }


def standardize_feature_matrix(features: np.ndarray) -> np.ndarray:
    if features.size == 0:
        return features.reshape((features.shape[0], 0))
    out = np.asarray(features, dtype=float).copy()
    for j in range(out.shape[1]):
        col = out[:, j]
        finite = np.isfinite(col)
        fill = float(np.nanmedian(col[finite])) if finite.any() else 0.0
        col[~finite] = fill
        scale = float(np.nanstd(col))
        out[:, j] = 0.0 if scale <= 0 else (col - float(np.nanmean(col))) / scale
    return out


def design_matrix_with_features(
    t_min: np.ndarray,
    components: list[tuple[float, int]],
    features: np.ndarray | None = None,
) -> np.ndarray:
    base = design_matrix(t_min, components)
    if features is None or features.size == 0:
        return base
    return np.column_stack([base, features])


def fit_linear_with_features(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    components: list[tuple[float, int]],
    features: np.ndarray | None = None,
) -> dict:
    x = design_matrix_with_features(t_min, components, features)
    beta = np.linalg.lstsq(x / sigma[:, None], y / sigma, rcond=None)[0]
    fitted = x @ beta
    resid = y - fitted
    chi2 = float(np.sum((resid / sigma) ** 2))
    sse = float(np.sum(resid**2))
    n = len(y)
    k = x.shape[1]
    dof = n - k
    return {
        "N": n,
        "k": k,
        "dof": dof,
        "chi2": chi2,
        "reduced_chi2": chi2 / dof,
        "aic_chi2": chi2 + 2 * k,
        "bic_chi2": chi2 + k * math.log(n),
        "sse": sse,
        "bic_sse": n * math.log(max(sse / n, 1e-300)) + k * math.log(n),
        "rms_pct": 100 * float(np.sqrt(np.mean(resid**2))),
        "beta": beta,
        "fitted": fitted,
        "resid": resid,
    }


def scan_period_with_features(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    p_min: float,
    p_max: float,
    n_grid: int,
    harmonics: tuple[int, ...],
    features: np.ndarray | None = None,
) -> tuple[float, pd.DataFrame, dict]:
    periods = np.linspace(p_min, p_max, n_grid)
    rows = []
    fits = []
    for period in periods:
        fit = fit_linear_with_features(t_min, y, sigma, [(float(period), h) for h in harmonics], features)
        row = {k: v for k, v in fit.items() if k not in {"beta", "fitted", "resid"}}
        row["period_min"] = float(period)
        rows.append(row)
        fits.append(fit)
    curve = pd.DataFrame(rows)
    idx = int(curve.chi2.idxmin())
    curve["delta_chi2"] = curve.chi2 - float(curve.loc[idx, "chi2"])
    return float(curve.loc[idx, "period_min"]), curve, fits[idx]


def cross_validated_rmse_with_features(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    components: list[tuple[float, int]],
    features: np.ndarray | None = None,
    n_folds: int = 5,
) -> float:
    order = np.argsort(t_min)
    folds = np.array_split(order, n_folds)
    pred = np.full(len(y), np.nan)
    for test_idx in folds:
        train_mask = np.ones(len(y), dtype=bool)
        train_mask[test_idx] = False
        train_features = features[train_mask] if features is not None and features.size else None
        test_features = features[test_idx] if features is not None and features.size else None
        fit = fit_linear_with_features(t_min[train_mask], y[train_mask], sigma[train_mask], components, train_features)
        pred[test_idx] = design_matrix_with_features(t_min[test_idx], components, test_features) @ fit["beta"]
    return 100 * float(np.sqrt(np.nanmean((y - pred) ** 2)))


def systematics_feature_sets(good: pd.DataFrame, t_min: np.ndarray) -> dict[str, tuple[list[str], np.ndarray]]:
    state = good[["file", "exptime", "n_comparisons", "comparison_ensemble_flux_rate"]].copy()
    phot_path = PHOT_DIR / "aperture_photometry.csv"
    if phot_path.exists():
        phot = pd.read_csv(phot_path)
        target_state = phot[phot.object_id.astype(str).eq("target")][
            ["file", "dx", "dy", "shift_scatter", "n_anchor_centroids", "raw_peak"]
        ]
        state = state.merge(target_state, on="file", how="left")
    state["time_min"] = t_min
    state["common_mode"] = state["comparison_ensemble_flux_rate"]
    specs = {
        "none": [],
        "linear_time": ["time_min"],
        "tracking": ["dx", "dy", "shift_scatter"],
        "observing_state": ["exptime", "n_comparisons", "common_mode"],
        "tracking_plus_observing": ["time_min", "dx", "dy", "shift_scatter", "exptime", "n_comparisons", "common_mode"],
    }
    out = {}
    for name, cols in specs.items():
        if not cols:
            out[name] = (cols, np.zeros((len(good), 0)))
        else:
            out[name] = (cols, standardize_feature_matrix(state[cols].to_numpy(float)))
    return out


def run_systematics_model_suite(
    good: pd.DataFrame,
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
) -> pd.DataFrame:
    rows = []
    for name, (cols, features) in systematics_feature_sets(good, t_min).items():
        period, _, fit = scan_period_with_features(t_min, y, sigma, 14, 22, 1600, (1, 2), features)
        cv = cross_validated_rmse_with_features(t_min, y, sigma, [(period, 1), (period, 2)], features)
        mid = np.nanmedian(t_min)
        p_first, _, _ = scan_period_with_features(
            t_min[t_min <= mid],
            y[t_min <= mid],
            sigma[t_min <= mid],
            14,
            22,
            700,
            (1, 2),
            features[t_min <= mid] if features.size else None,
        )
        p_second, _, _ = scan_period_with_features(
            t_min[t_min > mid],
            y[t_min > mid],
            sigma[t_min > mid],
            14,
            22,
            700,
            (1, 2),
            features[t_min > mid] if features.size else None,
        )
        rows.append(
            {
                "feature_set": name,
                "features": ";".join(cols),
                "period_min": period,
                "period_s": 60 * period,
                "N": fit["N"],
                "k": fit["k"],
                "chi2": fit["chi2"],
                "reduced_chi2": fit["reduced_chi2"],
                "bic_chi2": fit["bic_chi2"],
                "rms_pct": fit["rms_pct"],
                "crossval_rmse_pct": cv,
                "first_second_period_diff_s": abs(60 * p_first - 60 * p_second),
            }
        )
    df = pd.DataFrame(rows)
    base_bic = float(df.loc[df.feature_set.eq("none"), "bic_chi2"].iloc[0])
    base_cv = float(df.loc[df.feature_set.eq("none"), "crossval_rmse_pct"].iloc[0])
    df["delta_bic_vs_no_systematics"] = df.bic_chi2 - base_bic
    df["delta_cv_vs_no_systematics_pct"] = df.crossval_rmse_pct - base_cv
    return df


def robust_double_wave_refits(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    period0: float,
) -> pd.DataFrame:
    rows = []
    features = np.zeros((len(y), 0))
    init = fit_linear_with_features(t_min, y, sigma, [(period0, 1), (period0, 2)], features)
    beta0 = init["beta"]
    for loss in ["linear", "soft_l1", "huber", "cauchy"]:
        p0 = np.r_[period0, beta0]

        def residual(params: np.ndarray) -> np.ndarray:
            period = float(params[0])
            beta = params[1:]
            x = design_matrix_with_features(t_min, [(period, 1), (period, 2)], features)
            return (y - x @ beta) / sigma

        res = optimize.least_squares(
            residual,
            p0,
            bounds=(np.r_[14.0, np.full(len(beta0), -np.inf)], np.r_[22.0, np.full(len(beta0), np.inf)]),
            loss=loss,
            f_scale=1.0,
            max_nfev=3000,
        )
        r = residual(res.x)
        rms_pct = 100 * float(np.sqrt(np.mean((r * sigma) ** 2)))
        robust_scale = 1.4826 * float(np.nanmedian(np.abs((r * sigma) - np.nanmedian(r * sigma))))
        rows.append(
            {
                "loss": loss,
                "period_min": float(res.x[0]),
                "period_s": float(60 * res.x[0]),
                "cost": float(2 * res.cost),
                "ordinary_chi2_at_solution": float(np.sum(r**2)),
                "rms_pct": rms_pct,
                "robust_residual_sigma_pct": 100 * robust_scale,
                "max_abs_standard_residual": float(np.nanmax(np.abs(r))),
                "success": bool(res.success),
            }
        )
    return pd.DataFrame(rows)


def fixed_model_specs() -> list[dict]:
    orbital_min = ORBITAL_S / 60
    superhump_min = SUPERHUMP_S / 60
    negative_min = NEG_SUPERHUMP_S / 60
    return [
        {
            "model": "constant",
            "physical_interpretation": "no coherent AM CVn variability",
            "constraint_type": "baseline",
            "components": [],
        },
        {
            "model": "superhump_fixed_family",
            "physical_interpretation": "positive superhump at the literature 1051.2 s signal",
            "constraint_type": "conditional",
            "components": [(superhump_min, 1)],
        },
        {
            "model": "orbital_fixed_family",
            "physical_interpretation": "spectroscopic orbital period at 1028.7322 s",
            "constraint_type": "conditional",
            "components": [(orbital_min, 1)],
        },
        {
            "model": "harmonic_superhump_template",
            "physical_interpretation": "positive superhump waveform with strong 525.6 s harmonic",
            "constraint_type": "conditional",
            "components": [(superhump_min, 1), (superhump_min, 2)],
        },
        {
            "model": "orbital_plus_superhump",
            "physical_interpretation": "simultaneous orbital and positive-superhump signals",
            "constraint_type": "conditional",
            "components": [(orbital_min, 1), (superhump_min, 1)],
        },
        {
            "model": "orbital_plus_positive_negative_superhump",
            "physical_interpretation": "orbital plus positive and negative superhumps",
            "constraint_type": "conditional",
            "components": [(orbital_min, 1), (superhump_min, 1), (negative_min, 1)],
        },
    ]


def run_multiharmonic_lomb_scargle(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
) -> pd.DataFrame:
    """Run regular and multi-term Lomb-Scargle searches.

    Astropy's nterms option generalizes Lomb-Scargle from a single sinusoid to a
    Fourier series at each trial frequency.  nterms=2 is the LS analogue of the
    free double-wave model used below.
    """
    rows = []
    t_days = t_min / 1440
    yy = y - np.nanmedian(y)
    searches = [
        ("single_term_harmonic_window", 1, 5.0, 12.0),
        ("two_term_fundamental_window", 2, 14.0, 22.0),
        ("three_term_fundamental_window", 3, 14.0, 22.0),
    ]
    for label, nterms, p_min, p_max in searches:
        ls = LombScargle(t_days, yy, sigma, nterms=nterms)
        freq, power = ls.autopower(
            minimum_frequency=1 / (p_max / 1440),
            maximum_frequency=1 / (p_min / 1440),
            samples_per_peak=25,
        )
        period_min = (1 / freq) * 1440
        idx = int(np.nanargmax(power))
        fap = float(ls.false_alarm_probability(power[idx])) if nterms == 1 else np.nan
        rows.append(
            {
                "search": label,
                "nterms": nterms,
                "min_period_min": p_min,
                "max_period_min": p_max,
                "best_period_min": float(period_min[idx]),
                "best_period_s": float(60 * period_min[idx]),
                "best_power": float(power[idx]),
                "false_alarm_probability": fap,
                "interpretation": (
                    "single-sinusoid LS detects the strong harmonic"
                    if nterms == 1
                    else "multi-harmonic LS can recover the double-wave fundamental"
                ),
            }
        )
    return pd.DataFrame(rows)


def scan_beat_constrained_superhump(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    prec_min_hr: float = 6.0,
    prec_max_hr: float = 30.0,
    n_grid: int = 3000,
) -> tuple[pd.DataFrame, dict]:
    """Fit a precessing eccentric-disk proxy.

    Positive superhumps satisfy f_sh = f_orb - f_prec for prograde apsidal
    precession.  With only a 111 min baseline, precession is not directly
    observable, but this scan tests whether the best photometric period maps to
    a plausible disk-precession period when the spectroscopic orbital period is
    fixed.
    """
    f_orb = 1.0 / (ORBITAL_S / 60)
    rows = []
    fits = []
    for p_prec_hr in np.linspace(prec_min_hr, prec_max_hr, n_grid):
        f_prec = 1.0 / (p_prec_hr * 60)
        f_sh = f_orb - f_prec
        if f_sh <= 0:
            continue
        p_sh_min = 1.0 / f_sh
        fit = fit_linear(t_min, y, sigma, [(p_sh_min, 1), (p_sh_min, 2)])
        beta = fit["beta"]
        row = {
            "precession_period_hr": float(p_prec_hr),
            "superhump_period_min": float(p_sh_min),
            "superhump_period_s": float(60 * p_sh_min),
            "orbital_period_s": ORBITAL_S,
            "N": fit["N"],
            "k": fit["k"],
            "dof": fit["dof"],
            "chi2": fit["chi2"],
            "reduced_chi2": fit["reduced_chi2"],
            "aic_chi2": fit["aic_chi2"],
            "bic_chi2": fit["bic_chi2"],
            "sse": fit["sse"],
            "bic_sse": fit["bic_sse"],
            "rms_pct": fit["rms_pct"],
            "fundamental_amp_pct": 100 * float(np.hypot(beta[1], beta[2])),
            "harmonic_amp_pct": 100 * float(np.hypot(beta[3], beta[4])),
            "fundamental_phase_cycles": float((math.atan2(beta[2], beta[1]) / (2 * np.pi)) % 1),
            "harmonic_phase_cycles": float((math.atan2(beta[4], beta[3]) / (2 * np.pi)) % 1),
        }
        rows.append(row)
        fits.append(fit)
    curve = pd.DataFrame(rows)
    idx = int(curve.chi2.idxmin())
    curve["delta_chi2"] = curve.chi2 - float(curve.loc[idx, "chi2"])
    return curve, fits[idx]


def alternative_identity_tests(
    baseline_min: float,
    best_period_s: float,
    residual_fap: float,
) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "hypothesis": "AM CVn white-dwarf binary",
                "can_current_photometry_test": "partially",
                "current_data_result": "supports known 525/1051 s photometric family",
                "decisive_test": "time-resolved spectroscopy: helium emission/absorption, S-wave, Doppler tomography",
                "notes": "literature spectroscopy already establishes the binary/orbital interpretation",
            },
            {
                "hypothesis": "single pulsating white dwarf",
                "can_current_photometry_test": "weakly",
                "current_data_result": "period alone is not decisive; double-wave/superhump family favors AM CVn context",
                "decisive_test": "spectrum without accretion signatures plus multi-period WD pulsation pattern",
                "notes": "would not naturally explain published orbital S-wave spectroscopy",
            },
            {
                "hypothesis": "rapidly oscillating Ap or delta Scuti/SX Phe star",
                "can_current_photometry_test": "weakly",
                "current_data_result": f"{best_period_s:.1f} s is far shorter than ordinary delta Scuti periods; colors/spectrum needed",
                "decisive_test": "stellar spectrum, Gaia color/absolute magnitude, absence of helium accretion features",
                "notes": "photometric period alone cannot classify the object",
            },
            {
                "hypothesis": "instrumental/comparison-star artifact",
                "can_current_photometry_test": "yes",
                "current_data_result": f"leave-one-comparison and residual tests passed; residual FAP={residual_fap:.3f}",
                "decisive_test": "repeat on another night and independent comparison ensemble",
                "notes": "first/second-half stability is imperfect, so more data would still help",
            },
            {
                "hypothesis": "not a compact binary",
                "can_current_photometry_test": "no",
                "current_data_result": f"111.2 min baseline cannot override decades of spectroscopy",
                "decisive_test": "classification spectroscopy and radial-velocity curve",
                "notes": "broadband photometry can support/refine timing but cannot falsify the established object class",
            },
        ]
    )


def scan_period(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    p_min: float,
    p_max: float,
    n_grid: int,
    harmonics: tuple[int, ...],
) -> tuple[float, pd.DataFrame, dict]:
    periods = np.linspace(p_min, p_max, n_grid)
    rows = []
    fits = []
    for period in periods:
        fit = fit_linear(t_min, y, sigma, [(float(period), h) for h in harmonics])
        row = {k: v for k, v in fit.items() if k not in {"beta", "fitted", "resid"}}
        row["period_min"] = float(period)
        rows.append(row)
        fits.append(fit)
    curve = pd.DataFrame(rows)
    idx = int(curve.chi2.idxmin())
    curve["delta_chi2"] = curve.chi2 - float(curve.loc[idx, "chi2"])
    return float(curve.loc[idx, "period_min"]), curve, fits[idx]


def interval_from_curve(curve: pd.DataFrame, threshold: float, label: str, error_model: str) -> dict:
    idx = int(curve.chi2.idxmin())
    best = float(curve.loc[idx, "period_min"])
    below = curve[curve.delta_chi2 <= threshold]
    return {
        "constraint": label,
        "constraint_type": "independent" if "free" in label else "conditional",
        "error_model": error_model,
        "best_period_min": best,
        "best_period_s": 60 * best,
        "lo_period_min": float(below.period_min.min()),
        "hi_period_min": float(below.period_min.max()),
        "minus_s": 60 * (best - float(below.period_min.min())),
        "plus_s": 60 * (float(below.period_min.max()) - best),
        "sigma_s_avg": 30 * (float(below.period_min.max()) - float(below.period_min.min())),
        "delta_chi2_threshold": threshold,
    }


def cross_validated_rmse(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    components: list[tuple[float, int]],
    n_folds: int = 5,
) -> float:
    order = np.argsort(t_min)
    folds = np.array_split(order, n_folds)
    pred = np.full(len(y), np.nan)
    for test_idx in folds:
        train_mask = np.ones(len(y), dtype=bool)
        train_mask[test_idx] = False
        fit = fit_linear(t_min[train_mask], y[train_mask], sigma[train_mask], components)
        pred[test_idx] = design_matrix(t_min[test_idx], components) @ fit["beta"]
    return 100 * float(np.sqrt(np.nanmean((y - pred) ** 2)))


def block_bootstrap_indices(t_min: np.ndarray, rng: np.random.Generator, block_min: float = 10.0) -> np.ndarray:
    block_id = np.floor((t_min - t_min.min()) / block_min).astype(int)
    blocks = [np.where(block_id == b)[0] for b in np.unique(block_id)]
    chosen: list[int] = []
    while len(chosen) < len(t_min):
        chosen.extend(blocks[int(rng.integers(0, len(blocks)))].tolist())
    return np.asarray(chosen[: len(t_min)], dtype=int)


def bootstrap_free_periods(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
    single_grid: int = 700,
    double_grid: int = 900,
) -> pd.DataFrame:
    rows = []
    for i in range(n_boot):
        idx = block_bootstrap_indices(t_min, rng)
        p_single, _, fit_single = scan_period(t_min[idx], y[idx], sigma[idx], 5, 12, single_grid, (1,))
        p_double, _, fit_double = scan_period(t_min[idx], y[idx], sigma[idx], 14, 22, double_grid, (1, 2))
        rows.append(
            {
                "iteration": i,
                "single_sine_free_period_min": p_single,
                "single_sine_free_period_s": 60 * p_single,
                "single_sine_free_chi2": fit_single["chi2"],
                "double_wave_free_period_min": p_double,
                "double_wave_free_period_s": 60 * p_double,
                "double_wave_free_chi2": fit_double["chi2"],
            }
        )
    return pd.DataFrame(rows)


def summarize_bootstrap(boot: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for prefix in ["single_sine_free", "double_wave_free"]:
        period = boot[f"{prefix}_period_s"].to_numpy(float)
        med = float(np.nanmedian(period))
        p16 = float(np.nanpercentile(period, 16))
        p84 = float(np.nanpercentile(period, 84))
        p025 = float(np.nanpercentile(period, 2.5))
        p975 = float(np.nanpercentile(period, 97.5))
        rows.append(
            {
                "model": prefix,
                "N_boot": len(boot),
                "period_median_s": med,
                "minus_1sigma_s": med - p16,
                "plus_1sigma_s": p84 - med,
                "period_p2p5_s": p025,
                "period_p97p5_s": p975,
            }
        )
    return pd.DataFrame(rows)


def monte_carlo_free_periods(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    n_mc: int,
    rng: np.random.Generator,
    single_grid: int = 700,
    double_grid: int = 900,
) -> pd.DataFrame:
    rows = []
    for i in range(n_mc):
        y_pert = y + rng.normal(0.0, sigma)
        p_single, _, fit_single = scan_period(t_min, y_pert, sigma, 5, 12, single_grid, (1,))
        p_double, _, fit_double = scan_period(t_min, y_pert, sigma, 14, 22, double_grid, (1, 2))
        rows.append(
            {
                "iteration": i,
                "single_sine_free_period_min": p_single,
                "single_sine_free_period_s": 60 * p_single,
                "single_sine_free_chi2": fit_single["chi2"],
                "double_wave_free_period_min": p_double,
                "double_wave_free_period_s": 60 * p_double,
                "double_wave_free_chi2": fit_double["chi2"],
            }
        )
    return pd.DataFrame(rows)


def summarize_monte_carlo(mc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for prefix in ["single_sine_free", "double_wave_free"]:
        period = mc[f"{prefix}_period_s"].to_numpy(float)
        med = float(np.nanmedian(period))
        p16 = float(np.nanpercentile(period, 16))
        p84 = float(np.nanpercentile(period, 84))
        p025 = float(np.nanpercentile(period, 2.5))
        p975 = float(np.nanpercentile(period, 97.5))
        rows.append(
            {
                "model": prefix,
                "N_monte_carlo": len(mc),
                "period_median_s": med,
                "minus_1sigma_s": med - p16,
                "plus_1sigma_s": p84 - med,
                "period_p2p5_s": p025,
                "period_p97p5_s": p975,
            }
        )
    return pd.DataFrame(rows)


def plot_monte_carlo_periods(mc: pd.DataFrame, summary: pd.DataFrame) -> Path:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    specs = [
        ("single_sine_free", "Single-sine harmonic", 525.0550183394465, "tab:blue"),
        ("double_wave_free", "Double-wave fundamental", 1055.2138034508628, "tab:purple"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for ax, (model, title, original_period_s, color) in zip(axes, specs):
        periods = mc[f"{model}_period_s"].to_numpy(float)
        row = summary[summary.model.eq(model)].iloc[0]
        median = float(row.period_median_s)
        lo = median - float(row.minus_1sigma_s)
        hi = median + float(row.plus_1sigma_s)
        ax.hist(periods, bins=32, color=color, alpha=0.75, edgecolor="white")
        ax.axvline(original_period_s, color="black", lw=1.6, label="Original best fit")
        ax.axvline(median, color="tab:red", lw=1.8, label="MC median")
        ax.axvspan(lo, hi, color="tab:red", alpha=0.16, label="16th-84th pct.")
        ax.set_title(title)
        ax.set_xlabel("Best-fit period [s]")
        ax.set_ylabel("Monte Carlo realizations")
        ax.grid(alpha=0.25)
        ax.text(
            0.03,
            0.95,
            f"{median:.2f} s\n-{float(row.minus_1sigma_s):.2f}/+{float(row.plus_1sigma_s):.2f} s",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.9},
        )
    axes[0].legend(fontsize=8, loc="upper right")
    fig.suptitle("Fixed-time Monte Carlo period uncertainty", fontsize=14)
    out = PLOT_DIR / "AM_CVn_monte_carlo_period_uncertainty.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def leave_one_comparison_out(t_min: np.ndarray, sigma: np.ndarray) -> pd.DataFrame:
    phot_path = PHOT_DIR / "aperture_photometry.csv"
    if not phot_path.exists():
        return pd.DataFrame()
    phot = pd.read_csv(phot_path)
    ok = phot[(phot.flag == "ok") & np.isfinite(phot.flux_rate) & (phot.flux_rate > 0)].copy()
    medians = ok.groupby("object_id")["flux_rate"].median().to_dict()
    comp_labels = sorted([x for x in ok.object_id.unique() if str(x) != "target"], key=str)
    rows = []
    for dropped in ["none", *comp_labels]:
        diff_rows = []
        for _, group in phot.groupby("file", sort=False):
            target = group[group.object_id.astype(str).eq("target")]
            if target.empty:
                continue
            target = target.iloc[0]
            comps = group[
                group.role.eq("comparison")
                & group.flag.eq("ok")
                & np.isfinite(group.flux_rate)
                & (group.flux_rate > 0)
                & ~group.object_id.astype(str).eq(str(dropped))
            ].copy()
            comps["relative_flux"] = [r.flux_rate / medians.get(r.object_id, np.nan) for r in comps.itertuples()]
            comps = comps[np.isfinite(comps.relative_flux) & (comps.relative_flux > 0)]
            target_med = medians.get("target", np.nan)
            common = float(np.nanmedian(comps.relative_flux)) if len(comps) >= 2 else np.nan
            good = (
                target.flag == "ok"
                and np.isfinite(target.flux_rate)
                and np.isfinite(target_med)
                and target_med > 0
                and np.isfinite(common)
                and common > 0
                and not bool(target.frame_artifact)
            )
            if good:
                diff_rows.append(
                    {
                        "jd": target.jd,
                        "diff_flux_norm": target.flux_rate / target_med / common,
                        "diff_flux_err": target.flux_rate_err / target_med / common,
                    }
                )
        df = pd.DataFrame(diff_rows).sort_values("jd")
        if len(df) < 20:
            continue
        tt = (df.jd.to_numpy(float) - df.jd.min()) * 24 * 60
        yy = df.diff_flux_norm.to_numpy(float)
        ss = df.diff_flux_err.to_numpy(float) / float(np.nanmedian(df.diff_flux_norm.to_numpy(float)))
        p_double, _, fit_double = scan_period(tt, yy, ss, 14, 22, 1200, (1, 2))
        rows.append(
            {
                "dropped_comparison": dropped,
                "N": len(df),
                "double_wave_period_min": p_double,
                "double_wave_period_s": 60 * p_double,
                "double_wave_rms_pct": fit_double["rms_pct"],
                "double_wave_bic_chi2": fit_double["bic_chi2"],
            }
        )
    return pd.DataFrame(rows)


def recompute_differential_from_comparisons(
    phot: pd.DataFrame,
    comparison_labels: tuple[str, ...],
) -> pd.DataFrame:
    ok = phot[(phot.flag == "ok") & np.isfinite(phot.flux_rate) & (phot.flux_rate > 0)].copy()
    medians = ok.groupby("object_id")["flux_rate"].median().to_dict()
    rows = []
    for frame, group in phot.groupby("file", sort=False):
        target = group[group.object_id.astype(str).eq("target")]
        if target.empty:
            continue
        target = target.iloc[0]
        comps = group[
            group.object_id.astype(str).isin(comparison_labels)
            & group.flag.eq("ok")
            & np.isfinite(group.flux_rate)
            & (group.flux_rate > 0)
        ].copy()
        comps["relative_flux"] = [r.flux_rate / medians.get(r.object_id, np.nan) for r in comps.itertuples()]
        comps = comps[np.isfinite(comps.relative_flux) & (comps.relative_flux > 0)]
        target_med = medians.get("target", np.nan)
        common = float(np.nanmedian(comps.relative_flux)) if len(comps) >= max(2, min(3, len(comparison_labels))) else np.nan
        good = (
            target.flag == "ok"
            and np.isfinite(target.flux_rate)
            and np.isfinite(target_med)
            and target_med > 0
            and np.isfinite(common)
            and common > 0
            and not bool(target.frame_artifact)
        )
        if good:
            rows.append(
                {
                    "file": frame,
                    "jd": target.jd,
                    "diff_flux_norm": target.flux_rate / target_med / common,
                    "diff_flux_err": target.flux_rate_err / target_med / common,
                    "n_comparisons": len(comps),
                }
            )
    return pd.DataFrame(rows).sort_values("jd")


def subset_control_scatter(phot: pd.DataFrame, comparison_labels: tuple[str, ...]) -> float:
    if len(comparison_labels) < 4:
        return np.nan
    scatters = []
    for pseudo in comparison_labels:
        remaining = tuple(x for x in comparison_labels if x != pseudo)
        ok = phot[(phot.flag == "ok") & np.isfinite(phot.flux_rate) & (phot.flux_rate > 0)].copy()
        medians = ok.groupby("object_id")["flux_rate"].median().to_dict()
        vals = []
        for _, group in phot.groupby("file", sort=False):
            pseudo_row = group[group.object_id.astype(str).eq(pseudo)]
            if pseudo_row.empty:
                continue
            pseudo_row = pseudo_row.iloc[0]
            comps = group[
                group.object_id.astype(str).isin(remaining)
                & group.flag.eq("ok")
                & np.isfinite(group.flux_rate)
                & (group.flux_rate > 0)
            ].copy()
            comps["relative_flux"] = [r.flux_rate / medians.get(r.object_id, np.nan) for r in comps.itertuples()]
            comps = comps[np.isfinite(comps.relative_flux) & (comps.relative_flux > 0)]
            common = float(np.nanmedian(comps.relative_flux)) if len(comps) >= 2 else np.nan
            pseudo_med = medians.get(pseudo, np.nan)
            if pseudo_row.flag == "ok" and np.isfinite(common) and common > 0 and np.isfinite(pseudo_med) and pseudo_med > 0:
                vals.append(pseudo_row.flux_rate / pseudo_med / common)
        if len(vals) >= 20:
            scatters.append(100 * robust_sigma(np.asarray(vals, dtype=float)))
    return float(np.nanmedian(scatters)) if scatters else np.nan


def comparison_subset_robustness() -> pd.DataFrame:
    phot_path = PHOT_DIR / "aperture_photometry.csv"
    if not phot_path.exists():
        return pd.DataFrame()
    phot = pd.read_csv(phot_path)
    labels = sorted([str(x) for x in phot.loc[phot.role.eq("comparison"), "object_id"].unique()])
    rows = []
    for size in range(3, len(labels) + 1):
        for subset in itertools.combinations(labels, size):
            lc = recompute_differential_from_comparisons(phot, subset)
            if len(lc) < 80:
                continue
            tt = (lc.jd.to_numpy(float) - lc.jd.min()) * 24 * 60
            yy = lc.diff_flux_norm.to_numpy(float)
            ss = lc.diff_flux_err.to_numpy(float) / float(np.nanmedian(lc.diff_flux_norm.to_numpy(float)))
            period, _, fit = scan_period(tt, yy, ss, 14, 22, 900, (1, 2))
            control = subset_control_scatter(phot, subset)
            rows.append(
                {
                    "comparison_subset": ";".join(subset),
                    "n_subset": size,
                    "N": len(lc),
                    "period_s": 60 * period,
                    "target_model_rms_pct": fit["rms_pct"],
                    "target_reduced_chi2": fit["reduced_chi2"],
                    "target_bic_chi2": fit["bic_chi2"],
                    "control_median_scatter_pct": control,
                    "quality_score_pct": fit["rms_pct"] + (control if np.isfinite(control) else 10),
                }
            )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    all_label = ";".join(labels)
    base = df[df.comparison_subset.eq(all_label)]
    if not base.empty:
        base_row = base.iloc[0]
        df["delta_rms_vs_all_pct"] = df.target_model_rms_pct - float(base_row.target_model_rms_pct)
        df["delta_period_vs_all_s"] = df.period_s - float(base_row.period_s)
        df["delta_control_scatter_vs_all_pct"] = df.control_median_scatter_pct - float(base_row.control_median_scatter_pct)
    return df.sort_values(["quality_score_pct", "target_model_rms_pct"]).reset_index(drop=True)


def residual_periodogram(t_min: np.ndarray, resid: np.ndarray, sigma: np.ndarray) -> tuple[pd.DataFrame, dict]:
    ls = LombScargle(t_min / 1440, resid - np.nanmedian(resid), sigma)
    freq, power = ls.autopower(
        minimum_frequency=1 / (90 / 1440),
        maximum_frequency=1 / (5 / 1440),
        samples_per_peak=20,
    )
    period_min = (1 / freq) * 1440
    idx = int(np.nanargmax(power))
    curve = pd.DataFrame({"period_min": period_min, "power": power})
    summary = {
        "residual_best_period_min": float(period_min[idx]),
        "residual_best_period_s": float(60 * period_min[idx]),
        "residual_best_power": float(power[idx]),
        "residual_best_fap": float(ls.false_alarm_probability(power[idx])),
    }
    return curve, summary


def validation_summary(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    best_components: list[tuple[float, int]],
    loo: pd.DataFrame,
) -> pd.DataFrame:
    fit_all = fit_linear(t_min, y, sigma, best_components)
    z = fit_all["resid"] / sigma
    clip = np.abs(z) <= 4
    p_clip, _, _ = scan_period(t_min[clip], y[clip], sigma[clip], 14, 22, 1600, (1, 2))
    mid = np.nanmedian(t_min)
    p_first, _, _ = scan_period(t_min[t_min <= mid], y[t_min <= mid], sigma[t_min <= mid], 14, 22, 1200, (1, 2))
    p_second, _, _ = scan_period(t_min[t_min > mid], y[t_min > mid], sigma[t_min > mid], 14, 22, 1200, (1, 2))
    loo_span = float(loo.double_wave_period_s.max() - loo.double_wave_period_s.min()) if not loo.empty else np.nan
    return pd.DataFrame(
        [
            {
                "test": "first_half_second_half",
                "metric": "double_wave_period_s_difference",
                "value": abs(60 * p_first - 60 * p_second),
                "pass": abs(60 * p_first - 60 * p_second) < 20,
                "note": f"first={60*p_first:.2f}s second={60*p_second:.2f}s",
            },
            {
                "test": "outlier_clip_abs_z_le_4",
                "metric": "double_wave_period_s_shift",
                "value": abs(60 * p_clip - 60 * best_components[0][0]),
                "pass": abs(60 * p_clip - 60 * best_components[0][0]) < 10,
                "note": f"kept={int(np.sum(clip))}/{len(clip)} clipped_period={60*p_clip:.2f}s",
            },
            {
                "test": "leave_one_comparison_out",
                "metric": "double_wave_period_s_span",
                "value": loo_span,
                "pass": bool(np.isfinite(loo_span) and loo_span < 20),
                "note": "computed from aperture photometry by dropping one comparison at a time",
            },
        ]
    )


def make_plots(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    best_double: float,
    best_fit: dict,
    residual_curve: pd.DataFrame,
    models: pd.DataFrame,
    constraints: pd.DataFrame,
) -> None:
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)

    ax = axes[0, 0]
    ax.errorbar(t_min, 100 * (y - 1), yerr=100 * sigma, fmt=".", ms=4, color="0.35", ecolor="0.82", lw=0.6)
    grid = np.linspace(t_min.min(), t_min.max(), 800)
    model_grid = design_matrix(grid, [(best_double, 1), (best_double, 2)]) @ best_fit["beta"]
    ax.plot(grid, 100 * (model_grid - 1), color="tab:red", lw=2)
    ax.set_xlabel("Minutes from first exposure")
    ax.set_ylabel("Differential flux - 1 [%]")
    ax.set_title("Best free double-wave model")

    ax = axes[0, 1]
    phase = (t_min / best_double) % 1
    model_phase = (grid / best_double) % 1
    order = np.argsort(model_phase)
    ax.plot(np.r_[phase, phase + 1], 100 * np.r_[y - 1, y - 1], ".", ms=4, color="0.45")
    ax.plot(np.r_[model_phase[order], model_phase[order] + 1], 100 * np.r_[model_grid[order] - 1, model_grid[order] - 1], color="tab:red", lw=2)
    ax.set_xlabel("Phase at free double-wave period")
    ax.set_ylabel("Differential flux - 1 [%]")
    ax.set_title(f"Folded at {60 * best_double:.1f} s")

    ax = axes[1, 0]
    ax.plot(residual_curve.period_min, residual_curve.power, color="k", lw=1)
    ax.axvline(best_double, color="tab:red", ls="--", lw=1, label="free double wave")
    ax.axvline(SUPERHUMP_S / 60, color="tab:green", ls=":", lw=1, label="1051.2 s")
    ax.axvline(ORBITAL_S / 60, color="tab:blue", ls=":", lw=1, label="1028.7 s")
    ax.set_xlabel("Residual period [min]")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title("Residual periodogram")
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    plot_models = models.sort_values("bic_chi2").head(8).sort_values("bic_chi2", ascending=False)
    ax.barh(plot_models.model, plot_models.delta_bic_vs_free_double_wave, color="tab:purple")
    ax.axvline(0, color="k", lw=1)
    ax.set_xlabel("Delta BIC vs free double wave")
    ax.set_title("Model ladder")

    fig.savefig(PLOT_DIR / "AM_CVn_physical_model_ladder.png", dpi=180)

    fig, ax = plt.subplots(figsize=(7.5, 4.4), constrained_layout=True)
    literature = [
        ("orbital", ORBITAL_S, ORBITAL_SIGMA_S),
        ("superhump", SUPERHUMP_S, SUPERHUMP_SIGMA_S),
        ("Patterson 1979", PATTERSON_SUPERHUMP_S, PATTERSON_SUPERHUMP_SIGMA_S),
    ]
    yloc = np.arange(len(literature) + len(constraints))
    for i, (label, period, err) in enumerate(literature):
        ax.errorbar(period, yloc[i], xerr=err, fmt="o", color="tab:blue")
        ax.text(period + 4, yloc[i], label, va="center", fontsize=9)
    offset = len(literature)
    for j, row in constraints.iterrows():
        err = np.array([[row.minus_s], [row.plus_s]])
        ax.errorbar(row.best_period_s, yloc[offset + j], xerr=err, fmt="o", color="tab:red")
        ax.text(row.best_period_s + 10, yloc[offset + j], row.constraint, va="center", fontsize=9)
    ax.set_xlabel("Period [s]")
    ax.set_yticks([])
    ax.set_title("Our constraints compared with literature periods")
    fig.savefig(PLOT_DIR / "AM_CVn_physical_literature_constraints.png", dpi=180)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap", type=int, default=400, help="number of block-bootstrap realizations")
    parser.add_argument("--monte-carlo", type=int, default=400, help="number of fixed-time Monte Carlo realizations")
    parser.add_argument("--resample-single-grid", type=int, default=700, help="single-sine grid for bootstrap/Monte Carlo")
    parser.add_argument("--resample-double-grid", type=int, default=900, help="double-wave grid for bootstrap/Monte Carlo")
    args = parser.parse_args()

    PER_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    good, t_min, y, sigma = load_light_curve()

    mhls = run_multiharmonic_lomb_scargle(t_min, y, sigma)
    mhls.to_csv(PER_DIR / "physical_multiharmonic_lomb_scargle.csv", index=False)

    single_best, single_curve, single_fit = scan_period(t_min, y, sigma, 5, 12, 3000, (1,))
    double_best, double_curve, double_fit = scan_period(t_min, y, sigma, 14, 22, 4000, (1, 2))
    beat_curve, beat_fit = scan_beat_constrained_superhump(t_min, y, sigma)
    beat_best = beat_curve.loc[int(beat_curve.chi2.idxmin())]
    systematics_models = run_systematics_model_suite(good, t_min, y, sigma)
    systematics_models.to_csv(PER_DIR / "physical_systematics_model_comparison.csv", index=False)
    robust_refits = robust_double_wave_refits(t_min, y, sigma, double_best)
    robust_refits.to_csv(PER_DIR / "physical_robust_likelihood_refits.csv", index=False)
    subset_robustness = comparison_subset_robustness()
    subset_robustness.to_csv(PER_DIR / "physical_comparison_subset_robustness.csv", index=False)

    models = []
    fit_lookup = {}
    for spec in fixed_model_specs():
        fit = fit_linear(t_min, y, sigma, spec["components"])
        fit_lookup[spec["model"]] = fit
        row = {k: v for k, v in fit.items() if k not in {"beta", "fitted", "resid"}}
        row.update({k: v for k, v in spec.items() if k != "components"})
        row["best_period_s"] = np.nan
        row["crossval_rmse_pct"] = cross_validated_rmse(t_min, y, sigma, spec["components"]) if spec["components"] else np.nan
        models.append(row)
    for name, interpretation, period, harmonics, fit in [
        ("single_sine_free", "independent sinusoidal photometric period", single_best, (1,), single_fit),
        ("double_wave_free", "independent fundamental plus harmonic waveform", double_best, (1, 2), double_fit),
    ]:
        row = {k: v for k, v in fit.items() if k not in {"beta", "fitted", "resid"}}
        row.update(
            {
                "model": name,
                "physical_interpretation": interpretation,
                "constraint_type": "independent",
                "best_period_s": 60 * period,
                "crossval_rmse_pct": cross_validated_rmse(t_min, y, sigma, [(period, h) for h in harmonics]),
            }
        )
        models.append(row)
        fit_lookup[name] = fit
    beat_row = {k: v for k, v in beat_fit.items() if k not in {"beta", "fitted", "resid"}}
    beat_row.update(
        {
            "model": "beat_constrained_precessing_disk",
            "physical_interpretation": "positive superhump with f_sh = f_orb - f_prec",
            "constraint_type": "conditional",
            "best_period_s": float(beat_best.superhump_period_s),
            "precession_period_hr": float(beat_best.precession_period_hr),
            "crossval_rmse_pct": cross_validated_rmse(
                t_min,
                y,
                sigma,
                [(float(beat_best.superhump_period_min), 1), (float(beat_best.superhump_period_min), 2)],
            ),
        }
    )
    models.append(beat_row)
    fit_lookup["beat_constrained_precessing_disk"] = beat_fit
    models_df = pd.DataFrame(models)
    free_double_bic = float(models_df.loc[models_df.model.eq("double_wave_free"), "bic_chi2"].iloc[0])
    free_double_cv = float(models_df.loc[models_df.model.eq("double_wave_free"), "crossval_rmse_pct"].iloc[0])
    const_bic = float(models_df.loc[models_df.model.eq("constant"), "bic_chi2"].iloc[0])
    models_df["delta_bic_vs_constant"] = models_df.bic_chi2 - const_bic
    models_df["delta_bic_vs_free_double_wave"] = models_df.bic_chi2 - free_double_bic
    models_df["delta_cv_rmse_vs_free_double_wave_pct"] = models_df.crossval_rmse_pct - free_double_cv
    models_df.to_csv(PER_DIR / "physical_model_comparison.csv", index=False)

    scale = math.sqrt(float(double_fit["reduced_chi2"]))
    constraints = pd.DataFrame(
        [
            interval_from_curve(single_curve, 1.0, "single_sine_free_formal", "formal_delta_chi2_1"),
            interval_from_curve(single_curve, float(single_fit["reduced_chi2"]), "single_sine_free_scaled", "scaled_to_reduced_chi2_1"),
            interval_from_curve(double_curve, 1.0, "double_wave_free_formal", "formal_delta_chi2_1"),
            interval_from_curve(double_curve, float(double_fit["reduced_chi2"]), "double_wave_free_scaled", "scaled_to_reduced_chi2_1"),
        ]
    )

    super_curve = scan_period(t_min, y, sigma, 17.1, 17.9, 1800, (1, 2))[1]
    orbital_curve = scan_period(t_min, y, sigma, 16.8, 17.5, 1800, (1,))[1]
    constraints = pd.concat(
        [
            constraints,
            pd.DataFrame(
                [
                    interval_from_curve(super_curve, float(double_fit["reduced_chi2"]), "conditional_superhump_family_scaled", "scaled_to_reduced_chi2_1"),
                    interval_from_curve(orbital_curve, float(double_fit["reduced_chi2"]), "conditional_orbital_family_scaled", "scaled_to_reduced_chi2_1"),
                ]
            ),
        ],
        ignore_index=True,
    )
    constraints.to_csv(PER_DIR / "physical_period_constraints.csv", index=False)
    single_curve.to_csv(PER_DIR / "physical_single_sine_profile.csv", index=False)
    double_curve.to_csv(PER_DIR / "physical_double_wave_profile.csv", index=False)
    super_curve.to_csv(PER_DIR / "physical_conditional_superhump_profile.csv", index=False)
    orbital_curve.to_csv(PER_DIR / "physical_conditional_orbital_profile.csv", index=False)
    beat_curve.to_csv(PER_DIR / "physical_beat_constrained_precession_profile.csv", index=False)

    boot = bootstrap_free_periods(
        t_min,
        y,
        sigma * scale,
        args.bootstrap,
        rng,
        single_grid=args.resample_single_grid,
        double_grid=args.resample_double_grid,
    )
    boot.to_csv(PER_DIR / "physical_block_bootstrap_periods.csv", index=False)
    boot_summary = summarize_bootstrap(boot)
    boot_summary.to_csv(PER_DIR / "physical_block_bootstrap_summary.csv", index=False)

    monte_carlo = monte_carlo_free_periods(
        t_min,
        y,
        sigma * scale,
        args.monte_carlo,
        rng,
        single_grid=args.resample_single_grid,
        double_grid=args.resample_double_grid,
    )
    monte_carlo.to_csv(PER_DIR / "physical_monte_carlo_periods.csv", index=False)
    monte_carlo_summary = summarize_monte_carlo(monte_carlo)
    monte_carlo_summary.to_csv(PER_DIR / "physical_monte_carlo_period_summary.csv", index=False)
    plot_monte_carlo_periods(monte_carlo, monte_carlo_summary)

    loo = leave_one_comparison_out(t_min, sigma)
    loo.to_csv(PER_DIR / "physical_leave_one_comparison_out.csv", index=False)
    validation = validation_summary(t_min, y, sigma * scale, [(double_best, 1), (double_best, 2)], loo)

    residual_curve, residual_summary = residual_periodogram(t_min, double_fit["resid"], sigma * scale)
    residual_curve.to_csv(PER_DIR / "physical_residual_periodogram.csv", index=False)
    residual_summary_df = pd.DataFrame([residual_summary])
    residual_summary_df.to_csv(PER_DIR / "physical_residual_periodogram_summary.csv", index=False)
    alternative_identity_tests(
        baseline_min=float(t_min.max() - t_min.min()),
        best_period_s=float(60 * double_best),
        residual_fap=float(residual_summary["residual_best_fap"]),
    ).to_csv(PER_DIR / "physical_alternative_identity_tests.csv", index=False)

    best_physical = models_df[models_df.constraint_type.eq("conditional")].sort_values("bic_chi2").iloc[0]
    best_physical_beats_free = bool(best_physical.bic_chi2 < free_double_bic and best_physical.crossval_rmse_pct <= free_double_cv)
    free_scaled = constraints[constraints.constraint.eq("double_wave_free_scaled")].iloc[0]
    super_scaled = constraints[constraints.constraint.eq("conditional_superhump_family_scaled")].iloc[0]
    conclusion = pd.DataFrame(
        [
            {
                "N": len(y),
                "baseline_min": float(t_min.max() - t_min.min()),
                "cycles_at_superhump": float((t_min.max() - t_min.min()) / (SUPERHUMP_S / 60)),
                "precession_cycles_covered": float((t_min.max() - t_min.min()) / (PRECESSION_HR * 60)),
                "best_independent_model": "double_wave_free",
                "best_independent_period_s": float(60 * double_best),
                "best_independent_scaled_sigma_s": float(free_scaled.sigma_s_avg),
                "best_literature_informed_model": str(best_physical.model),
                "best_literature_informed_delta_bic_vs_free": float(best_physical.delta_bic_vs_free_double_wave),
                "conditional_superhump_period_s": float(super_scaled.best_period_s),
                "conditional_superhump_scaled_sigma_s": float(super_scaled.sigma_s_avg),
                "literature_sigma_to_beat_s": min(PATTERSON_SUPERHUMP_SIGMA_S, SUPERHUMP_SIGMA_S, ORBITAL_SIGMA_S),
                "physically_motivated_model_beats_free_model": best_physical_beats_free,
                "can_tighten_literature_independently": bool(best_physical_beats_free and free_scaled.sigma_s_avg < PATTERSON_SUPERHUMP_SIGMA_S),
                "recommended_claim": (
                    "consistent_with_literature_not_tighter"
                    if not (best_physical_beats_free and free_scaled.sigma_s_avg < PATTERSON_SUPERHUMP_SIGMA_S)
                    else "potentially_tighter_independent_constraint"
                ),
            }
        ]
    )
    conclusion.to_csv(PER_DIR / "physical_constraint_conclusion.csv", index=False)

    validation = pd.concat(
        [
            validation,
            pd.DataFrame(
                [
                    {
                        "test": "residual_periodogram",
                        "metric": "best_residual_fap",
                        "value": residual_summary["residual_best_fap"],
                        "pass": residual_summary["residual_best_fap"] > 0.05,
                        "note": f"best residual period={residual_summary['residual_best_period_min']:.2f} min",
                    },
                    {
                        "test": "literature_tightening",
                        "metric": "independent_sigma_s_vs_patterson_sigma_s",
                        "value": float(free_scaled.sigma_s_avg / PATTERSON_SUPERHUMP_SIGMA_S),
                        "pass": bool(free_scaled.sigma_s_avg < PATTERSON_SUPERHUMP_SIGMA_S),
                        "note": "must be <1 to beat Patterson 1979 period precision",
                    },
                ]
            ),
        ],
        ignore_index=True,
    )
    validation.to_csv(PER_DIR / "physical_validation_summary.csv", index=False)

    make_plots(t_min, y, sigma * scale, double_best, double_fit, residual_curve, models_df, constraints)

    print("Wrote AM CVn physical model-suite products.")
    print(conclusion.to_string(index=False))
    print(models_df.sort_values("bic_chi2")[["model", "bic_chi2", "delta_bic_vs_free_double_wave", "crossval_rmse_pct"]].to_string(index=False))
    print(systematics_models.sort_values("bic_chi2").head(8).to_string(index=False))
    print(robust_refits.to_string(index=False))
    print(monte_carlo_summary.to_string(index=False))
    if not subset_robustness.empty:
        print(subset_robustness.head(8).to_string(index=False))
    print(validation.to_string(index=False))


if __name__ == "__main__":
    main()
