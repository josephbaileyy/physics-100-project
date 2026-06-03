"""Exploratory AM CVn period modeling and literature comparison.

This script keeps the measurement data-driven: literature periods are used only
after fitting our light curve. Outputs go to the existing AM CVn periodogram and
plot folders.
"""
from __future__ import annotations

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
from scipy import stats
from sklearn.mixture import GaussianMixture


ROOT = Path(__file__).resolve().parents[2]
PHOT_DIR = ROOT / "analysis/am_cvn/photometry/sequence_X42421BZ"
PER_DIR = ROOT / "analysis/am_cvn/periodograms/sequence_X42421BZ"
PLOT_DIR = ROOT / "analysis/am_cvn/plots"
N_BOOT = 10_000
RNG_SEED = 20260529
SAMPLES_PER_PEAK = 20
HALF_FAMILY_MAX_MIN = 12.0
LONG_FAMILY_MIN_MIN = 12.0
LONG_FAMILY_MAX_MIN = 25.0


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


def design_matrix(t_min: np.ndarray, period_min: float | None, model: str) -> np.ndarray:
    if model == "constant":
        return np.ones((len(t_min), 1))
    if period_min is None:
        raise ValueError("period_min is required")
    phase = 2 * np.pi * t_min / period_min
    if model == "sine":
        return np.column_stack([np.ones(len(t_min)), np.sin(phase), np.cos(phase)])
    if model == "double_wave":
        return np.column_stack(
            [np.ones(len(t_min)), np.sin(phase), np.cos(phase), np.sin(2 * phase), np.cos(2 * phase)]
        )
    raise ValueError(f"Unknown model: {model}")


def weighted_fit(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    period_min: float | None,
    model: str,
) -> dict:
    x = design_matrix(t_min, period_min, model)
    beta = np.linalg.lstsq(x / sigma[:, None], y / sigma, rcond=None)[0]
    fitted = x @ beta
    resid = y - fitted
    chi2 = float(np.sum((resid / sigma) ** 2))
    sse = float(np.sum(resid**2))
    n = len(y)
    k = x.shape[1]
    dof = n - k
    row = {
        "period_min": np.nan if period_min is None else period_min,
        "model_type": model,
        "N": n,
        "k": k,
        "dof": dof,
        "chi2": chi2,
        "reduced_chi2": chi2 / dof,
        "chi2_p_value": float(stats.chi2.sf(chi2, dof)),
        "bic_chi2": chi2 + k * math.log(n),
        "sse": sse,
        "bic_sse": n * math.log(max(sse / n, 1e-300)) + k * math.log(n),
        "rms_pct": 100 * float(np.sqrt(np.mean(resid**2))),
    }
    if k >= 3:
        row["fundamental_amp_pct"] = 100 * float(np.hypot(beta[1], beta[2]))
    if k >= 5:
        row["harmonic_amp_pct"] = 100 * float(np.hypot(beta[3], beta[4]))
    return row | {"beta": beta, "fitted": fitted, "resid": resid}


def fit_without_arrays(*args, **kwargs) -> dict:
    row = weighted_fit(*args, **kwargs)
    row.pop("beta")
    row.pop("fitted")
    row.pop("resid")
    return row


def chi2_curve(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    model: str,
    p_min: float,
    p_max: float,
    n_grid: int = 6000,
) -> pd.DataFrame:
    periods = np.linspace(p_min, p_max, n_grid)
    rows = [fit_without_arrays(t_min, y, sigma, float(p), model) for p in periods]
    out = pd.DataFrame(rows)
    out["delta_chi2"] = out.chi2 - float(out.chi2.min())
    return out


def interval_from_curve(curve: pd.DataFrame, threshold: float, label: str, error_model: str) -> dict:
    idx = int(curve.chi2.idxmin())
    best = float(curve.loc[idx, "period_min"])
    below = curve[curve.delta_chi2 <= threshold]
    return {
        "experiment": label,
        "error_model": error_model,
        "best_period_min": best,
        "best_period_s": 60 * best,
        "min_chi2": float(curve.loc[idx, "chi2"]),
        "delta_chi2_threshold": threshold,
        "lo_period_min": float(below.period_min.min()),
        "hi_period_min": float(below.period_min.max()),
        "minus_min": best - float(below.period_min.min()),
        "plus_min": float(below.period_min.max()) - best,
        "minus_s": 60 * (best - float(below.period_min.min())),
        "plus_s": 60 * (float(below.period_min.max()) - best),
    }


def run_ls_window(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    min_period: float,
    max_period: float,
    label: str,
) -> tuple[dict, pd.DataFrame]:
    t_days = t_min / (24 * 60)
    ls = LombScargle(t_days, y - np.nanmedian(y), sigma)
    freq, power = ls.autopower(
        minimum_frequency=1 / (max_period / (24 * 60)),
        maximum_frequency=1 / (min_period / (24 * 60)),
        samples_per_peak=SAMPLES_PER_PEAK,
    )
    period = (1 / freq) * 24 * 60
    idx = int(np.nanargmax(power))
    best_period = float(period[idx])
    best_power = float(power[idx])
    halfmax = best_power / 2
    above = np.where(power >= halfmax)[0]
    local = above[(above >= max(0, idx - 500)) & (above <= min(len(power) - 1, idx + 500))]
    width = float(np.nanmax(period[local]) - np.nanmin(period[local])) if len(local) else np.nan
    row = {
        "window": label,
        "min_period_min": min_period,
        "max_period_min": max_period,
        "best_period_min": best_period,
        "best_period_s": 60 * best_period,
        "best_power": best_power,
        "false_alarm_probability": float(ls.false_alarm_probability(best_power)),
        "half_power_width_min": width,
        "family": classify_period(best_period),
    }
    curve = pd.DataFrame({"period_min": period, "frequency_cpd": freq, "power": power, "window": label})
    return row, curve


def classify_period(period_min: float) -> str:
    if period_min < 0.8:
        return "very_fast_26s_candidate"
    if 1.7 <= period_min <= 2.2:
        return "115_120s_candidate"
    if period_min < HALF_FAMILY_MAX_MIN:
        return "half_period_family"
    if LONG_FAMILY_MIN_MIN <= period_min <= LONG_FAMILY_MAX_MIN:
        return "long_period_family"
    return "other"


def f_test(parent: dict, child: dict) -> dict:
    df_num = child["k"] - parent["k"]
    df_den = child["dof"]
    f_stat = ((parent["chi2"] - child["chi2"]) / df_num) / (child["chi2"] / df_den)
    return {
        "parent_model": parent["model"],
        "child_model": child["model"],
        "df_num": df_num,
        "df_den": df_den,
        "F": float(f_stat),
        "p_value": float(stats.f.sf(f_stat, df_num, df_den)),
        "delta_chi2": parent["chi2"] - child["chi2"],
        "delta_bic_chi2": child["bic_chi2"] - parent["bic_chi2"],
        "delta_bic_sse": child["bic_sse"] - parent["bic_sse"],
    }


def summarize_bootstrap(df: pd.DataFrame, label: str) -> dict:
    periods = df.best_period_min.to_numpy(float)
    faps = df.false_alarm_probability.to_numpy(float)
    return {
        "experiment": label,
        "N": len(df),
        "period_median_min": float(np.nanmedian(periods)),
        "period_p16_min": float(np.nanpercentile(periods, 16)),
        "period_p84_min": float(np.nanpercentile(periods, 84)),
        "period_p2p5_min": float(np.nanpercentile(periods, 2.5)),
        "period_p97p5_min": float(np.nanpercentile(periods, 97.5)),
        "fap_median": float(np.nanmedian(faps)),
        "fap_p16": float(np.nanpercentile(faps, 16)),
        "fap_p84": float(np.nanpercentile(faps, 84)),
        "fap_below_0p05_frac": float(np.mean(faps < 0.05)),
        "half_family_frac": float(np.mean(periods < HALF_FAMILY_MAX_MIN)),
        "long_family_frac": float(np.mean((periods >= LONG_FAMILY_MIN_MIN) & (periods <= LONG_FAMILY_MAX_MIN))),
    }


def family_intervals(df: pd.DataFrame, label: str) -> list[dict]:
    rows = []
    families = [
        ("half_period_family", df.best_period_min < HALF_FAMILY_MAX_MIN),
        ("long_period_family", (df.best_period_min >= LONG_FAMILY_MIN_MIN) & (df.best_period_min <= LONG_FAMILY_MAX_MIN)),
    ]
    for family, mask in families:
        sub = df[mask]
        if sub.empty:
            continue
        periods = sub.best_period_min.to_numpy(float)
        faps = sub.false_alarm_probability.to_numpy(float)
        med = float(np.nanmedian(periods))
        p16 = float(np.nanpercentile(periods, 16))
        p84 = float(np.nanpercentile(periods, 84))
        rows.append(
            {
                "bootstrap": label,
                "reported_component": family,
                "weight": len(sub) / len(df),
                "period_median_min": med,
                "minus_1sigma_min": med - p16,
                "plus_1sigma_min": p84 - med,
                "period_median_s": 60 * med,
                "minus_1sigma_s": 60 * (med - p16),
                "plus_1sigma_s": 60 * (p84 - med),
                "fap_median": float(np.nanmedian(faps)),
                "fap_below_0p05_frac": float(np.mean(faps < 0.05)),
            }
        )
        if family == "half_period_family":
            rows.append(
                {
                    "bootstrap": label,
                    "reported_component": "half_period_doubled_to_fundamental",
                    "weight": len(sub) / len(df),
                    "period_median_min": 2 * med,
                    "minus_1sigma_min": 2 * (med - p16),
                    "plus_1sigma_min": 2 * (p84 - med),
                    "period_median_s": 120 * med,
                    "minus_1sigma_s": 120 * (med - p16),
                    "plus_1sigma_s": 120 * (p84 - med),
                    "fap_median": float(np.nanmedian(faps)),
                    "fap_below_0p05_frac": float(np.mean(faps < 0.05)),
                }
            )
    return rows


def gmm_summary(df: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    p = df.best_period_min.to_numpy(float)
    for space, x in [("period_min", p), ("frequency_cpd", 1440 / p)]:
        gm = GaussianMixture(n_components=2, covariance_type="full", n_init=50, random_state=RNG_SEED)
        gm.fit(x.reshape(-1, 1))
        means = gm.means_.ravel()
        sigmas = np.sqrt(gm.covariances_.reshape(2))
        weights = gm.weights_.ravel()
        order = np.argsort(means if space == "period_min" else -means)
        for family, idx in zip(["half_period_family", "long_period_family"], order):
            if space == "period_min":
                mean_min = means[idx]
                sigma_min = sigmas[idx]
            else:
                mean_min = 1440 / means[idx]
                sigma_min = 1440 * sigmas[idx] / means[idx] ** 2
            rows.append(
                {
                    "bootstrap": label,
                    "fit_space": space,
                    "family": family,
                    "weight": weights[idx],
                    "model_period_mean_min": mean_min,
                    "model_period_sigma_min": sigma_min,
                    "model_period_mean_s": 60 * mean_min,
                    "model_period_sigma_s": 60 * sigma_min,
                }
            )
    return pd.DataFrame(rows)


def block_bootstrap_indices(t_min: np.ndarray, rng: np.random.Generator, block_min: float = 10.0) -> np.ndarray:
    block_id = np.floor((t_min - t_min.min()) / block_min).astype(int)
    blocks = [np.where(block_id == b)[0] for b in np.unique(block_id)]
    chosen: list[int] = []
    while len(chosen) < len(t_min):
        block = blocks[int(rng.integers(0, len(blocks)))]
        chosen.extend(block.tolist())
    return np.asarray(chosen[: len(t_min)], dtype=int)


def run_bootstrap_variants(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    scaled_sigma: np.ndarray,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    variants: dict[str, list[dict]] = {
        "row_resample": [],
        "perturb_formal_errors": [],
        "perturb_scaled_errors": [],
        "block_resample_10min": [],
    }
    n = len(y)
    for i in range(N_BOOT):
        idx = rng.integers(0, n, n)
        variants["row_resample"].append(run_ls_window(t_min[idx], y[idx], sigma[idx], 5, 90, "bootstrap")[0] | {"iteration": i})
        y_formal = y + rng.normal(0, sigma)
        variants["perturb_formal_errors"].append(run_ls_window(t_min, y_formal, sigma, 5, 90, "bootstrap")[0] | {"iteration": i})
        y_scaled = y + rng.normal(0, scaled_sigma)
        variants["perturb_scaled_errors"].append(
            run_ls_window(t_min, y_scaled, scaled_sigma, 5, 90, "bootstrap")[0] | {"iteration": i}
        )
        bidx = block_bootstrap_indices(t_min, rng)
        variants["block_resample_10min"].append(
            run_ls_window(t_min[bidx], y[bidx], sigma[bidx], 5, 90, "bootstrap")[0] | {"iteration": i}
        )

    all_rows = []
    summary_rows = []
    family_rows = []
    gmm_rows = []
    for label, rows in variants.items():
        df = pd.DataFrame(rows)
        df["experiment"] = label
        df.to_csv(PER_DIR / f"exploratory_bootstrap_{label}.csv", index=False)
        all_rows.append(df)
        summary_rows.append(summarize_bootstrap(df, label))
        family_rows.extend(family_intervals(df, label))
        gmm_rows.append(gmm_summary(df, label))
    return (
        pd.DataFrame(summary_rows),
        pd.DataFrame(family_rows),
        pd.concat(gmm_rows, ignore_index=True),
    )


def subsample_sweep(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    rng: np.random.Generator,
    n_iter: int = 2500,
) -> pd.DataFrame:
    rows = []
    n = len(y)
    for frac in [0.4, 0.6, 0.8, 1.0]:
        size = max(5, int(round(frac * n)))
        for i in range(n_iter):
            idx = np.sort(rng.choice(n, size=size, replace=False))
            row = run_ls_window(t_min[idx], y[idx], sigma[idx], 5, 90, "subsample")[0]
            row.update({"iteration": i, "sample_fraction": frac, "sample_size": size})
            rows.append(row)
    return pd.DataFrame(rows)


def model_comparison(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    ls_half: float,
    double_best: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, dict]]:
    model_specs = [
        ("constant_flux", "constant", None),
        ("sine_free_half_period", "sine", ls_half),
        ("sine_doubled_half_as_fundamental", "sine", 2 * ls_half),
        ("sine_double_wave_best_period", "sine", double_best),
        ("double_wave_best_period", "double_wave", double_best),
        ("sine_literature_1028p732s_orbital", "sine", 1028.7322 / 60),
        ("sine_literature_1051p2s_superhump", "sine", 1051.2 / 60),
    ]
    rows = []
    fits = {}
    for name, model, period in model_specs:
        fit = weighted_fit(t_min, y, sigma, period, model)
        fits[name] = fit
        row = {k: v for k, v in fit.items() if k not in {"beta", "fitted", "resid"}}
        row["model"] = name
        rows.append(row)
    df = pd.DataFrame(rows)
    const = df.loc[df.model.eq("constant_flux")].iloc[0]
    df["delta_bic_chi2_vs_constant"] = df.bic_chi2 - const.bic_chi2
    df["delta_bic_sse_vs_constant"] = df.bic_sse - const.bic_sse
    lookup = {row["model"]: row for row in df.to_dict("records")}
    ftests = pd.DataFrame(
        [
            f_test(lookup["constant_flux"], lookup["sine_free_half_period"]),
            f_test(lookup["constant_flux"], lookup["double_wave_best_period"]),
            f_test(lookup["sine_double_wave_best_period"], lookup["double_wave_best_period"]),
        ]
    )
    return df, ftests, fits


def time_local_phase_test(t_min: np.ndarray, y: np.ndarray, sigma: np.ndarray, period_min: float) -> pd.DataFrame:
    rows = []
    mid = np.nanmedian(t_min)
    for label, mask in [("first_half", t_min <= mid), ("second_half", t_min > mid), ("all", np.ones_like(t_min, dtype=bool))]:
        fit = weighted_fit(t_min[mask], y[mask], sigma[mask], period_min, "sine")
        beta = fit["beta"]
        amp = float(np.hypot(beta[1], beta[2]))
        phase_rad = float(math.atan2(beta[2], beta[1]))
        rows.append(
            {
                "segment": label,
                "N": int(np.sum(mask)),
                "period_min": period_min,
                "amp_pct": 100 * amp,
                "phase_rad": phase_rad,
                "phase_cycles": (phase_rad / (2 * np.pi)) % 1,
                "chi2": fit["chi2"],
                "reduced_chi2": fit["reduced_chi2"],
                "rms_pct": fit["rms_pct"],
            }
        )
    return pd.DataFrame(rows)


def detectability_summary(good: pd.DataFrame, t_min: np.ndarray) -> pd.DataFrame:
    duration = float(t_min.max() - t_min.min())
    cadence_s = float(np.nanmedian(np.diff(np.sort(t_min))) * 60)
    rows = []
    for label, period_s in [
        ("26.3s oscillation", 26.3),
        ("115s oscillation", 115.0),
        ("120s oscillation", 120.0),
        ("525s harmonic", 525.0),
        ("1028.7322s orbital", 1028.7322),
        ("1051.2s superhump", 1051.2),
        ("13.38h precession", 13.38 * 3600),
    ]:
        exptime = good.exptime.to_numpy(float)
        attenuation = np.abs(np.sinc(exptime / period_s))
        rows.append(
            {
                "signal": label,
                "period_s": period_s,
                "period_min": period_s / 60,
                "cycles_in_dataset": duration / (period_s / 60),
                "median_samples_per_cycle": period_s / cadence_s,
                "median_exposure_smearing_attenuation": float(np.nanmedian(attenuation)),
                "min_exposure_smearing_attenuation": float(np.nanmin(attenuation)),
                "detectability_note": detectability_note(period_s, duration, cadence_s, attenuation),
            }
        )
    return pd.DataFrame(rows)


def detectability_note(period_s: float, duration_min: float, cadence_s: float, attenuation: np.ndarray) -> str:
    if period_s < np.nanmax(attenuation * 0 + 45):
        return "not reliable: period is shorter than or comparable to exposures and below cadence Nyquist"
    if period_s / cadence_s < 3:
        return "marginal: only about 2-3 samples per cycle"
    if duration_min / (period_s / 60) < 1:
        return "not measurable: less than one cycle in the dataset"
    return "measurable as a photometric component if amplitude is high enough"


def literature_table(our_periods: dict[str, tuple[float, float]]) -> pd.DataFrame:
    refs = [
        {
            "paper": "Patterson et al. 1979",
            "period_type": "photometric superhump-like 1051s family",
            "period_s": 1051.212,
            "sigma_s": 0.015,
            "method": "extrema timing and O-C ephemeris",
            "classification": "same physical signal",
            "testability": "consistent, but our uncertainty is too large for dP/dt confirmation",
        },
        {
            "paper": "Kruszewski & Semeniuk 1992",
            "period_type": "dominant 525s harmonic / 1051s family",
            "period_s": 1051.2,
            "sigma_s": np.nan,
            "method": "1962 photometry reanalysis; period not uniquely determined",
            "classification": "same physical signal",
            "testability": "consistent",
        },
        {
            "paper": "Kruszewski & Semeniuk 1992",
            "period_type": "separate 1023s mode",
            "period_s": 1023.02,
            "sigma_s": np.nan,
            "method": "1962 photometry reanalysis",
            "classification": "different or uncertain mode",
            "testability": "our superhump-family signal is inconsistent, but this is not the same signal",
        },
        {
            "paper": "Nelemans et al. 2001",
            "period_type": "spectroscopic orbital period",
            "period_s": 1028.7322,
            "sigma_s": 0.0003,
            "method": "spectroscopic S-wave/orbital ephemeris",
            "classification": "different signal",
            "testability": "not detected as the dominant photometric signal",
        },
        {
            "paper": "Nelemans et al. 2001",
            "period_type": "positive superhump",
            "period_s": 1051.2,
            "sigma_s": np.nan,
            "method": "photometric superhump interpretation",
            "classification": "same physical signal",
            "testability": "consistent",
        },
        {
            "paper": "Nelemans et al. 2001",
            "period_type": "disk precession",
            "period_s": 13.38 * 3600,
            "sigma_s": np.nan,
            "method": "beat between orbital and superhump periods",
            "classification": "not measurable here",
            "testability": "not testable: our baseline covers about 0.14 cycles",
        },
        {
            "paper": "Roelofs et al. 2006",
            "period_type": "spectroscopic orbital period",
            "period_s": 1028.7322,
            "sigma_s": 0.0003,
            "method": "phase-binned spectroscopy/Doppler tomography",
            "classification": "different signal",
            "testability": "not detected as the dominant photometric signal",
        },
    ]
    rows = []
    for ref in refs:
        for our_label, (period_s, sigma_s) in our_periods.items():
            combined = math.sqrt(sigma_s**2 + (0 if not np.isfinite(ref["sigma_s"]) else ref["sigma_s"] ** 2))
            rows.append(
                ref
                | {
                    "our_measurement": our_label,
                    "our_period_s": period_s,
                    "our_sigma_s": sigma_s,
                    "difference_our_minus_ref_s": period_s - ref["period_s"],
                    "raw_z_using_our_sigma": (period_s - ref["period_s"]) / combined if combined > 0 else np.nan,
                    "z_if_same_signal": (period_s - ref["period_s"]) / combined
                    if combined > 0 and ref["classification"] == "same physical signal"
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)


def plot_window_comparison(window_curves: pd.DataFrame, window_summary: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11, 5), constrained_layout=True)
    for label, group in window_curves.groupby("window"):
        ax.plot(group.period_min, group.power, lw=0.9, label=label)
    for row in window_summary.itertuples(index=False):
        ax.axvline(row.best_period_min, color="0.2", alpha=0.15, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("Period [min]")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title("AM CVn Lomb-Scargle window sensitivity")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, ncols=2)
    fig.savefig(PLOT_DIR / "AM_CVn_exploratory_periodogram_windows.png", dpi=180)
    plt.close(fig)


def plot_folds(t_min: np.ndarray, y: np.ndarray, sigma: np.ndarray, periods: dict[str, float]) -> None:
    y_pct = 100 * (y - 1)
    yerr = 100 * sigma
    fig, axes = plt.subplots(3, 2, figsize=(12, 11), constrained_layout=True)
    axes = axes.ravel()
    for ax, (label, period) in zip(axes, periods.items()):
        phase = (t_min / period) % 1
        order = np.argsort(phase)
        ax.errorbar(phase[order], y_pct[order], yerr=yerr[order], fmt=".", ms=4, alpha=0.6, ecolor="0.82")
        ax.errorbar(phase[order] + 1, y_pct[order], yerr=yerr[order], fmt=".", ms=4, alpha=0.6, ecolor="0.82")
        model = "double_wave" if "double-wave" in label else "sine"
        fit = weighted_fit(t_min, y, sigma, period, model)
        grid = np.linspace(0, 2, 1200)
        pred = design_matrix(grid * period, period, model) @ fit["beta"]
        ax.plot(grid, 100 * (pred - 1), color="tab:red", lw=1.8)
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xlim(0, 2)
        ax.set_xlabel("Phase")
        ax.set_ylabel("Residual [%]")
        ax.set_title(label)
        ax.grid(alpha=0.25)
    if len(periods) < len(axes):
        for ax in axes[len(periods) :]:
            ax.axis("off")
    fig.savefig(PLOT_DIR / "AM_CVn_exploratory_folded_models.png", dpi=180)
    plt.close(fig)


def plot_bootstrap_histograms(summary_files: list[Path], gmm_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    axes = axes.ravel()
    xgrid = np.linspace(5, 25, 1200)
    for ax, path in zip(axes, summary_files):
        label = path.stem.replace("exploratory_bootstrap_", "")
        df = pd.read_csv(path)
        ax.hist(df.best_period_min, bins=np.linspace(5, 25, 90), color="tab:blue", alpha=0.7, density=True)
        sub = gmm_df[(gmm_df.bootstrap == label) & (gmm_df.fit_space == "period_min")]
        for row in sub.itertuples(index=False):
            pdf = row.weight * stats.norm.pdf(xgrid, row.model_period_mean_min, row.model_period_sigma_min)
            ax.plot(xgrid, pdf, lw=1.5, label=row.family)
        ax.axvline(1051.2 / 60, color="tab:green", ls="--", lw=1, label="1051.2 s")
        ax.axvline(1028.7322 / 60, color="tab:red", ls=":", lw=1, label="1028.7 s")
        ax.set_xlim(5, 25)
        ax.set_xlabel("Best LS period [min]")
        ax.set_ylabel("density")
        ax.set_title(label)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.savefig(PLOT_DIR / "AM_CVn_exploratory_bootstrap_family_histograms.png", dpi=180)
    plt.close(fig)


def plot_residual_periodograms(curves: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 5), constrained_layout=True)
    for label, group in curves.groupby("window"):
        ax.plot(group.period_min, group.power, lw=1.0, label=label)
    ax.set_xscale("log")
    ax.set_xlabel("Residual period [min]")
    ax.set_ylabel("Lomb-Scargle power")
    ax.set_title("Residual periodograms after removing superhump-family double wave")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.savefig(PLOT_DIR / "AM_CVn_exploratory_residual_periodograms.png", dpi=180)
    plt.close(fig)


def main() -> None:
    PER_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    good, t_min, y, sigma = load_light_curve()

    cadence_s = float(np.nanmedian(np.diff(np.sort(t_min))) * 60)
    baseline_min = float(t_min.max() - t_min.min())

    ls_windows = [
        ("short_5_to_12min", 5.0, 12.0),
        ("focused_7_to_10min", 7.0, 10.0),
        ("long_14_to_22min", 14.0, 22.0),
        ("wide_3_to_180min", 3.0, 180.0),
        ("literature_525s_harmonic", 8.3, 9.2),
        ("literature_1028p7s_orbital", 16.8, 17.5),
        ("literature_1051s_superhump", 17.1, 17.9),
        ("fast_115_to_120s", 1.85, 2.05),
        ("fast_26p3s", 0.38, 0.50),
    ]
    window_rows = []
    curve_rows = []
    for label, pmin, pmax in ls_windows:
        row, curve = run_ls_window(t_min, y, sigma, pmin, pmax, label)
        window_rows.append(row)
        curve_rows.append(curve)
    window_summary = pd.DataFrame(window_rows)
    window_curves = pd.concat(curve_rows, ignore_index=True)
    window_summary.to_csv(PER_DIR / "exploratory_lomb_scargle_window_summary.csv", index=False)
    window_curves.to_csv(PER_DIR / "exploratory_lomb_scargle_window_curves.csv", index=False)

    ls_half = float(window_summary.loc[window_summary.window.eq("focused_7_to_10min"), "best_period_min"].iloc[0])
    chi2_sine = chi2_curve(t_min, y, sigma, "sine", 7.0, 10.0)
    chi2_double = chi2_curve(t_min, y, sigma, "double_wave", 14.0, 20.0)
    chi2_sine.to_csv(PER_DIR / "exploratory_sine_period_chi2_curve_7_to_10min.csv", index=False)
    chi2_double.to_csv(PER_DIR / "exploratory_double_wave_period_chi2_curve_14_to_20min.csv", index=False)
    double_best = float(chi2_double.loc[chi2_double.chi2.idxmin(), "period_min"])

    models, ftests, fit_objects = model_comparison(t_min, y, sigma, ls_half, double_best)
    best_model = models.sort_values("bic_chi2").iloc[0]
    error_scale = float(np.sqrt(best_model.reduced_chi2))
    scaled_sigma = sigma * error_scale
    models.to_csv(PER_DIR / "exploratory_model_comparison.csv", index=False)
    ftests.to_csv(PER_DIR / "exploratory_model_f_tests.csv", index=False)

    intervals = []
    for label, curve, model_name in [
        ("sine_half_period_family", chi2_sine, "sine_free_half_period"),
        ("double_wave_long_period_family", chi2_double, "double_wave_best_period"),
    ]:
        red = float(models.loc[models.model.eq(model_name), "reduced_chi2"].iloc[0])
        intervals.append(interval_from_curve(curve, 1.0, label, "formal_delta_chi2_1"))
        intervals.append(interval_from_curve(curve, red, label, "scaled_to_reduced_chi2_1"))
    intervals_df = pd.DataFrame(intervals)
    intervals_df.to_csv(PER_DIR / "exploratory_period_delta_chi2_intervals.csv", index=False)

    boot_summary, boot_families, gmm = run_bootstrap_variants(t_min, y, sigma, scaled_sigma, rng)
    boot_summary.to_csv(PER_DIR / "exploratory_bootstrap_period_summary.csv", index=False)
    boot_families.to_csv(PER_DIR / "exploratory_bootstrap_family_intervals.csv", index=False)
    gmm.to_csv(PER_DIR / "exploratory_bootstrap_gmm_summary.csv", index=False)

    sweep = subsample_sweep(t_min, y, sigma, rng)
    sweep.to_csv(PER_DIR / "exploratory_subsample_sweep_lomb_scargle.csv", index=False)
    sweep_summary = (
        sweep.groupby("sample_fraction")
        .agg(
            N=("best_period_min", "size"),
            period_median_min=("best_period_min", "median"),
            period_p16_min=("best_period_min", lambda s: np.nanpercentile(s, 16)),
            period_p84_min=("best_period_min", lambda s: np.nanpercentile(s, 84)),
            fap_median=("false_alarm_probability", "median"),
            fap_below_0p05_frac=("false_alarm_probability", lambda s: np.mean(s < 0.05)),
            half_family_frac=("best_period_min", lambda s: np.mean(s < HALF_FAMILY_MAX_MIN)),
            long_family_frac=("best_period_min", lambda s: np.mean((s >= LONG_FAMILY_MIN_MIN) & (s <= LONG_FAMILY_MAX_MIN))),
        )
        .reset_index()
    )
    sweep_summary.to_csv(PER_DIR / "exploratory_subsample_sweep_summary.csv", index=False)

    resid = fit_objects["double_wave_best_period"]["resid"]
    residual_windows = [
        ("residual_5_to_90min", 5.0, 90.0),
        ("residual_115_to_120s", 1.85, 2.05),
        ("residual_26p3s", 0.38, 0.50),
    ]
    residual_rows = []
    residual_curves = []
    for label, pmin, pmax in residual_windows:
        row, curve = run_ls_window(t_min, resid + 1, sigma, pmin, pmax, label)
        residual_rows.append(row)
        residual_curves.append(curve)
    residual_summary = pd.DataFrame(residual_rows)
    residual_curve_df = pd.concat(residual_curves, ignore_index=True)
    residual_summary.to_csv(PER_DIR / "exploratory_residual_periodogram_summary.csv", index=False)
    residual_curve_df.to_csv(PER_DIR / "exploratory_residual_periodogram_curves.csv", index=False)

    local_phase = pd.concat(
        [
            time_local_phase_test(t_min, y, sigma, ls_half).assign(fit_family="half_period"),
            time_local_phase_test(t_min, y, sigma, double_best).assign(fit_family="long_period"),
        ],
        ignore_index=True,
    )
    local_phase.to_csv(PER_DIR / "exploratory_time_local_amplitude_phase.csv", index=False)

    detectability = detectability_summary(good, t_min)
    detectability.to_csv(PER_DIR / "exploratory_signal_detectability_summary.csv", index=False)

    our_periods = {
        "doubled_half_period_scaled_delta_chi2": (
            2 * float(intervals_df.loc[intervals_df.experiment.eq("sine_half_period_family") & intervals_df.error_model.eq("scaled_to_reduced_chi2_1"), "best_period_s"].iloc[0]),
            2
            * float(
                0.5
                * (
                    intervals_df.loc[
                        intervals_df.experiment.eq("sine_half_period_family")
                        & intervals_df.error_model.eq("scaled_to_reduced_chi2_1"),
                        "minus_s",
                    ].iloc[0]
                    + intervals_df.loc[
                        intervals_df.experiment.eq("sine_half_period_family")
                        & intervals_df.error_model.eq("scaled_to_reduced_chi2_1"),
                        "plus_s",
                    ].iloc[0]
                )
            ),
        ),
        "double_wave_scaled_delta_chi2": (
            float(
                intervals_df.loc[
                    intervals_df.experiment.eq("double_wave_long_period_family")
                    & intervals_df.error_model.eq("scaled_to_reduced_chi2_1"),
                    "best_period_s",
                ].iloc[0]
            ),
            float(
                0.5
                * (
                    intervals_df.loc[
                        intervals_df.experiment.eq("double_wave_long_period_family")
                        & intervals_df.error_model.eq("scaled_to_reduced_chi2_1"),
                        "minus_s",
                    ].iloc[0]
                    + intervals_df.loc[
                        intervals_df.experiment.eq("double_wave_long_period_family")
                        & intervals_df.error_model.eq("scaled_to_reduced_chi2_1"),
                        "plus_s",
                    ].iloc[0]
                )
            ),
        ),
    }
    literature = literature_table(our_periods)
    literature.to_csv(PER_DIR / "exploratory_literature_comparison.csv", index=False)

    diagnostics = pd.DataFrame(
        [
            {
                "N": len(y),
                "baseline_min": baseline_min,
                "median_cadence_s": cadence_s,
                "median_formal_error_pct": 100 * float(np.nanmedian(sigma)),
                "robust_scatter_pct": 100 * robust_sigma(y),
                "best_bic_model": str(best_model.model),
                "best_bic_reduced_chi2": float(best_model.reduced_chi2),
                "error_scale_from_best_bic_model": error_scale,
                "precession_cycles_covered_13p38h": baseline_min / (13.38 * 60),
            }
        ]
    )
    diagnostics.to_csv(PER_DIR / "exploratory_dataset_diagnostics.csv", index=False)

    plot_window_comparison(window_curves, window_summary)
    plot_folds(
        t_min,
        y,
        sigma,
        {
            f"half-period fit: {ls_half:.3f} min": ls_half,
            f"doubled half-period: {2 * ls_half:.3f} min": 2 * ls_half,
            f"double-wave best: {double_best:.3f} min": double_best,
            "orbital literature: 1028.7322 s": 1028.7322 / 60,
            "superhump literature: 1051.2 s": 1051.2 / 60,
        },
    )
    bootstrap_files = [
        PER_DIR / "exploratory_bootstrap_row_resample.csv",
        PER_DIR / "exploratory_bootstrap_perturb_formal_errors.csv",
        PER_DIR / "exploratory_bootstrap_perturb_scaled_errors.csv",
        PER_DIR / "exploratory_bootstrap_block_resample_10min.csv",
    ]
    plot_bootstrap_histograms(bootstrap_files, gmm)
    plot_residual_periodograms(residual_curve_df)

    print("Wrote exploratory AM CVn period modeling products.")
    print(diagnostics.to_string(index=False))
    print(models.sort_values("bic_chi2").head(8).to_string(index=False))
    print(boot_summary.to_string(index=False))
    print(detectability.to_string(index=False))


if __name__ == "__main__":
    main()
