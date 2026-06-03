"""AM CVn model diagnostics and period-uncertainty experiments."""
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


ROOT = Path(__file__).resolve().parents[2]
PHOT_DIR = ROOT / "analysis/am_cvn/photometry/sequence_X42421BZ"
PER_DIR = ROOT / "analysis/am_cvn/periodograms/sequence_X42421BZ"
PLOT_DIR = ROOT / "analysis/am_cvn/plots"
LS_MIN_PERIOD = 5.0
LS_MAX_PERIOD = 90.0
LITERATURE_PERIOD = 17.14
N_BOOT = 1000
RNG_SEED = 20260529


def robust_sigma(x: np.ndarray) -> float:
    med = np.nanmedian(x)
    return float(1.4826 * np.nanmedian(np.abs(x - med)))


def load_light_curve() -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray]:
    diff = pd.read_csv(PHOT_DIR / "differential_light_curve.csv")
    good = diff[diff.good.astype(bool) & np.isfinite(diff.diff_flux_norm)].sort_values("jd").copy()
    median_diff_flux = float(np.nanmedian(good.diff_flux.to_numpy(float)))
    t_min = (good.jd.to_numpy(float) - good.jd.min()) * 24 * 60
    y = good.diff_flux_norm.to_numpy(float)
    sigma = good.diff_flux_err.to_numpy(float) / median_diff_flux
    ok = np.isfinite(t_min) & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    return good.iloc[np.where(ok)[0]].copy(), t_min[ok], y[ok], sigma[ok]


def design_matrix(t_min: np.ndarray, period_min: float | None, model: str) -> np.ndarray:
    if model == "constant":
        return np.ones((len(t_min), 1))
    if period_min is None:
        raise ValueError("period_min is required for periodic models")
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
    out = {
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
        out["fundamental_amp_pct"] = 100 * float(np.hypot(beta[1], beta[2]))
    if k >= 5:
        out["harmonic_amp_pct"] = 100 * float(np.hypot(beta[3], beta[4]))
    return out


def chi2_curve(
    t_min: np.ndarray,
    y: np.ndarray,
    sigma: np.ndarray,
    model: str,
    p_min: float,
    p_max: float,
    n_grid: int = 5000,
) -> pd.DataFrame:
    periods = np.linspace(p_min, p_max, n_grid)
    rows = [weighted_fit(t_min, y, sigma, float(period), model) for period in periods]
    out = pd.DataFrame(rows)
    best = float(out.chi2.min())
    out["delta_chi2"] = out.chi2 - best
    return out


def threshold_interval(curve: pd.DataFrame, delta: float) -> dict:
    idx = int(curve.chi2.idxmin())
    best_period = float(curve.loc[idx, "period_min"])
    best_chi2 = float(curve.loc[idx, "chi2"])
    below = curve[curve.delta_chi2 <= delta]
    return {
        "best_period_min": best_period,
        "min_chi2": best_chi2,
        "delta_chi2_threshold": delta,
        "lo_period_min": float(below.period_min.min()),
        "hi_period_min": float(below.period_min.max()),
        "minus_min": best_period - float(below.period_min.min()),
        "plus_min": float(below.period_min.max()) - best_period,
    }


def run_lomb_scargle(t_min: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> dict:
    t_days = t_min / (24 * 60)
    ls = LombScargle(t_days, y - np.nanmedian(y), sigma)
    freq, power = ls.autopower(
        minimum_frequency=1 / (LS_MAX_PERIOD / (24 * 60)),
        maximum_frequency=1 / (LS_MIN_PERIOD / (24 * 60)),
        samples_per_peak=20,
    )
    idx = int(np.nanargmax(power))
    period_min = float((1 / freq[idx]) * 24 * 60)
    return {
        "best_period_min": period_min,
        "best_power": float(power[idx]),
        "false_alarm_probability": float(ls.false_alarm_probability(power[idx])),
    }


def summarize_bootstrap(df: pd.DataFrame, label: str) -> dict:
    periods = df.best_period_min.to_numpy(float)
    faps = df.false_alarm_probability.to_numpy(float)
    near_ls = np.abs(periods - np.nanmedian(periods)) < 1.0
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
        "period_near_median_pm1min_frac": float(np.mean(near_ls)),
    }


def main() -> None:
    PER_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RNG_SEED)
    good, t_min, y, sigma = load_light_curve()

    period_results = pd.read_csv(PER_DIR / "period_results.csv")
    ls_period = float(period_results.loc[period_results.run.eq("short_5_to_90min"), "best_period_min"].iloc[0])

    model_specs = [
        ("constant_flux", "constant", None),
        ("sine_LS_best", "sine", ls_period),
        ("sine_17p14min", "sine", LITERATURE_PERIOD),
        ("double_wave_17p14min", "double_wave", LITERATURE_PERIOD),
    ]
    fits = []
    for name, model, period in model_specs:
        row = weighted_fit(t_min, y, sigma, period, model)
        row["model"] = name
        fits.append(row)
    fits_df = pd.DataFrame(fits)
    const = fits_df.loc[fits_df.model.eq("constant_flux")].iloc[0]
    fits_df["delta_bic_chi2_vs_constant"] = fits_df.bic_chi2 - const.bic_chi2
    fits_df["delta_bic_sse_vs_constant"] = fits_df.bic_sse - const.bic_sse

    best_double = fits_df.loc[fits_df.model.eq("double_wave_17p14min")].iloc[0]
    error_scale = float(np.sqrt(best_double.reduced_chi2))

    diagnostics = {
        "N_good": len(good),
        "median_formal_error_pct": 100 * float(np.nanmedian(sigma)),
        "mean_formal_error_pct": 100 * float(np.nanmean(sigma)),
        "robust_target_scatter_pct": 100 * robust_sigma(y),
        "best_model": "double_wave_17p14min",
        "best_model_reduced_chi2": float(best_double.reduced_chi2),
        "suggested_error_scale_from_best_model": error_scale,
        "scaled_median_error_pct": 100 * float(np.nanmedian(sigma)) * error_scale,
    }

    chi2_sine = chi2_curve(t_min, y, sigma, "sine", 6.0, 12.0)
    chi2_double = chi2_curve(t_min, y, sigma, "double_wave", 14.0, 20.0)
    chi2_sine.to_csv(PER_DIR / "sine_period_chi2_curve_6_to_12min.csv", index=False)
    chi2_double.to_csv(PER_DIR / "double_wave_period_chi2_curve_14_to_20min.csv", index=False)

    intervals = []
    for label, curve, best_row in [
        ("sine_around_LS_peak", chi2_sine, fits_df.loc[fits_df.model.eq("sine_LS_best")].iloc[0]),
        ("double_wave_around_17min", chi2_double, best_double),
    ]:
        formal = threshold_interval(curve, 1.0)
        formal["experiment"] = label
        formal["error_model"] = "formal_delta_chi2_1"
        intervals.append(formal)
        scaled = threshold_interval(curve, float(best_row.reduced_chi2))
        scaled["experiment"] = label
        scaled["error_model"] = "scaled_to_reduced_chi2_1"
        intervals.append(scaled)
    interval_df = pd.DataFrame(intervals)
    interval_df.to_csv(PER_DIR / "period_delta_chi2_intervals.csv", index=False)

    boot_rows = []
    for i in range(N_BOOT):
        sample = rng.integers(0, len(y), len(y))
        row = run_lomb_scargle(t_min[sample], y[sample], sigma[sample])
        row["iteration"] = i
        boot_rows.append(row)
    boot_subset = pd.DataFrame(boot_rows)
    boot_subset.to_csv(PER_DIR / "bootstrap_resample_lomb_scargle.csv", index=False)

    perturb_rows = []
    for i in range(N_BOOT):
        y_pert = y + rng.normal(0.0, sigma * error_scale)
        row = run_lomb_scargle(t_min, y_pert, sigma * error_scale)
        row["iteration"] = i
        perturb_rows.append(row)
    boot_perturb = pd.DataFrame(perturb_rows)
    boot_perturb.to_csv(PER_DIR / "bootstrap_perturbed_lomb_scargle.csv", index=False)

    boot_summary = pd.DataFrame(
        [
            summarize_bootstrap(boot_subset, "row_resample_with_replacement"),
            summarize_bootstrap(boot_perturb, "perturb_by_scaled_errors"),
        ]
    )
    boot_summary.to_csv(PER_DIR / "bootstrap_period_summary.csv", index=False)

    pd.DataFrame([diagnostics]).to_csv(PER_DIR / "photometry_error_diagnostics.csv", index=False)
    fits_df.to_csv(PER_DIR / "model_chi2_f_bic_summary_refined.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    axes[0, 0].plot(chi2_sine.period_min, chi2_sine.delta_chi2, color="tab:red")
    axes[0, 0].axhline(1, color="0.4", ls="--", lw=1, label="Delta chi2 = 1")
    axes[0, 0].axhline(float(fits_df.loc[fits_df.model.eq("sine_LS_best"), "reduced_chi2"].iloc[0]), color="0.2", ls=":", lw=1, label="scaled threshold")
    axes[0, 0].set_xlabel("Sine period [min]")
    axes[0, 0].set_ylabel("Delta chi2")
    axes[0, 0].set_title("Sine period uncertainty")
    axes[0, 0].grid(alpha=0.25)
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(chi2_double.period_min, chi2_double.delta_chi2, color="tab:purple")
    axes[0, 1].axhline(1, color="0.4", ls="--", lw=1)
    axes[0, 1].axhline(float(best_double.reduced_chi2), color="0.2", ls=":", lw=1)
    axes[0, 1].set_xlabel("Double-wave fundamental period [min]")
    axes[0, 1].set_ylabel("Delta chi2")
    axes[0, 1].set_title("Double-wave period uncertainty")
    axes[0, 1].grid(alpha=0.25)

    axes[1, 0].hist(boot_subset.best_period_min, bins=40, color="tab:blue", alpha=0.8)
    axes[1, 0].axvline(ls_period, color="tab:red", lw=1.5, label="original LS")
    axes[1, 0].set_xlabel("Bootstrap LS best period [min]")
    axes[1, 0].set_ylabel("count")
    axes[1, 0].set_title("Row-resample bootstrap")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].hist(boot_perturb.best_period_min, bins=40, color="tab:green", alpha=0.8)
    axes[1, 1].axvline(ls_period, color="tab:red", lw=1.5, label="original LS")
    axes[1, 1].set_xlabel("Perturbed-data LS best period [min]")
    axes[1, 1].set_ylabel("count")
    axes[1, 1].set_title("Perturb by scaled errors")
    axes[1, 1].legend(fontsize=8)
    fig.savefig(PLOT_DIR / "AM_CVn_period_uncertainty_diagnostics.png", dpi=180)
    plt.close(fig)

    print("Diagnostics")
    for key, value in diagnostics.items():
        print(f"{key}: {value}")
    print("\nModel fits")
    print(fits_df.to_string(index=False))
    print("\nDelta-chi2 intervals")
    print(interval_df.to_string(index=False))
    print("\nBootstrap summary")
    print(boot_summary.to_string(index=False))


if __name__ == "__main__":
    main()
